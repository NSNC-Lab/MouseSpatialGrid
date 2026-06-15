"""
Benchmark E-prop vs PyTorch backprop for one cell.

The benchmark uses cell 7 parameters from LongRunResults.mat by default and
compares:
  - full PyTorch backprop through time on CPU and CUDA
  - truncated PyTorch backprop through time on CPU and CUDA
  - the single-cell generated E-prop solver on CPU
  - a CuPy version generated from the same single-cell E-prop solver

Outputs:
  - benchmark_results.json
  - memory_summary.csv
  - speed_summary.csv
  - benchmark_report.md
  - speed_fit_outputs.png
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import importlib.util
import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
except ModuleNotFoundError as exc:
    raise SystemExit("This benchmark needs PyTorch. Activate the SCC eprop/torch environment first.") from exc

try:
    import cupy as cp
except ModuleNotFoundError as exc:
    cp = None

try:
    from scipy.io import loadmat
except ModuleNotFoundError as exc:
    raise SystemExit("This benchmark needs scipy to load .mat files.") from exc

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None

try:
    import psutil
except ModuleNotFoundError:
    psutil = None

from run_single_cell_backprop import (
    EPROP_DIR,
    InitialParams,
    clean_argv,
    load_cell_raster,
    positive,
    surrogate_spike,
)


PARAM_NAMES = [
    "strf_gain",
    "strf_latency",
    "ron_g_inc",
    "on_ron_gsyn",
    "off_ron_gsyn",
    "on_sonoff_gsyn",
    "off_sonoff_gsyn",
    "sonoff_ron_gsyn",
]


@dataclass
class PhaseProfile:
    zero_grad_s: float = 0.0
    forward_s: float = 0.0
    backward_s: float = 0.0
    optimizer_s: float = 0.0
    best_trace_s: float = 0.0
    total_s: float = 0.0


@dataclass
class BenchResult:
    test: str
    method: str
    batch_size: int
    epochs: int
    final_loss: float
    best_loss: float
    best_batch_id: int
    best_epoch: int
    wall_time_s: float
    peak_cuda_allocated_mb: Optional[float]
    peak_cuda_reserved_mb: Optional[float]
    cupy_pool_used_mb: Optional[float]
    cupy_pool_reserved_mb: Optional[float]
    rss_start_mb: Optional[float]
    rss_end_mb: Optional[float]
    rss_delta_mb: Optional[float]
    phase_profile: Dict[str, float]
    memory_culprits: List[Dict[str, object]]
    phase_memory: Dict[str, float]


@contextlib.contextmanager
def eprop_import_path() -> Iterable[None]:
    sys.path.insert(0, str(EPROP_DIR))
    try:
        yield
    finally:
        with contextlib.suppress(ValueError):
            sys.path.remove(str(EPROP_DIR))


def sync_torch(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.synchronize(device)


def reset_torch_memory(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats(device)


def torch_memory_mb(device: torch.device) -> Tuple[Optional[float], Optional[float]]:
    if device.type != "cuda":
        return None, None
    sync_torch(device)
    return (
        torch.cuda.max_memory_allocated(device) / 1024**2,
        torch.cuda.max_memory_reserved(device) / 1024**2,
    )


def reset_cupy_memory() -> None:
    cp.cuda.Stream.null.synchronize()
    cp.get_default_memory_pool().free_all_blocks()
    cp.get_default_pinned_memory_pool().free_all_blocks()


def cupy_pool_mb() -> Tuple[float, float]:
    cp.cuda.Stream.null.synchronize()
    pool = cp.get_default_memory_pool()
    return pool.used_bytes() / 1024**2, pool.total_bytes() / 1024**2


def current_rss_mb() -> Optional[float]:
    if psutil is None:
        return None
    return psutil.Process(os.getpid()).memory_info().rss / 1024**2


def mb(num_bytes: float) -> float:
    return float(num_bytes) / 1024**2


def tensor_bytes(*arrays: torch.Tensor) -> int:
    return int(sum(a.numel() * a.element_size() for a in arrays))


def np_bytes(*arrays: np.ndarray) -> int:
    return int(sum(a.nbytes for a in arrays))


def inv_softplus_np(values: np.ndarray) -> np.ndarray:
    values = np.maximum(values.astype(np.float32), 1e-10)
    return np.log(np.expm1(values)).astype(np.float32)


def inv_sigmoid_range_np(values: np.ndarray, low: float, high: float) -> np.ndarray:
    values = np.asarray(values, dtype=np.float32)
    unit = (values - low) / (high - low)
    unit = np.clip(unit, 1e-6, 1.0 - 1e-6)
    return np.log(unit / (1.0 - unit)).astype(np.float32)


def sigmoid_range(raw: torch.Tensor, low: float, high: float) -> torch.Tensor:
    return low + (high - low) * torch.sigmoid(raw)


def param_vec_to_initial(param_vec: np.ndarray) -> InitialParams:
    param_vec = np.asarray(param_vec, dtype=np.float32).reshape(-1)
    return InitialParams(
        strf_gain=float(param_vec[0]),
        strf_latency=float(param_vec[1]),
        ron_g_inc=float(param_vec[2]),
        on_ron_gsyn=float(param_vec[3]),
        off_ron_gsyn=float(param_vec[4]),
        on_sonoff_gsyn=float(param_vec[5]),
        off_sonoff_gsyn=float(param_vec[6]),
        sonoff_ron_gsyn=float(param_vec[7]),
    )


def load_cell_params(long_run_mat: Path, cell_id: int, batch_size: int, perturb: bool = False) -> np.ndarray:
    mat = loadmat(long_run_mat, squeeze_me=True, struct_as_record=False)
    if "params" not in mat:
        raise KeyError(f"{long_run_mat} does not contain a 'params' variable.")

    params = np.asarray(mat["params"], dtype=np.float32)
    if params.ndim == 3:
        cell_params = params[:, cell_id - 1, :]
    elif params.ndim == 2:
        cell_params = params
    else:
        raise ValueError(f"Unexpected params shape: {params.shape}")

    repeats = int(math.ceil(batch_size / cell_params.shape[1]))
    out = np.tile(cell_params, repeats)[:, :batch_size].copy()
    if perturb:
        perturb_factors = np.ones((out.shape[0], 1), dtype=np.float32)
        perturb_factors[3] = 1.15
        perturb_factors[5] = 0.85
        perturb_factors[7] = 1.10
        out *= perturb_factors
    return out


def make_default_long_run_path() -> Path:
    return EPROP_DIR.parents[1] / "JointMDS_Space" / "LongRunResults.mat"


def generate_inputs(
    params_2d: np.ndarray,
    batch_size: int,
    spontaneous_fr: float,
    sim_len: int,
    seed: int,
) -> Dict[str, np.ndarray]:
    np.random.seed(seed)
    with eprop_import_path():
        from input_handler import call_inputs
        from strf_handler import call_strfs

        with clean_argv():
            target_dict = call_strfs(params_2d, batch_size, 1)
            spks = call_inputs(1, spontaneous_fr, batch_size, target_dict)

    on_key = "locs_masker_None_target_0_on"
    off_key = "locs_masker_None_target_0_off"
    on_raw = spks[on_key]["stimulus_0_poisson_spks"]
    off_raw = spks[off_key]["stimulus_0_poisson_spks"]
    noise_raw = spks["noise_masker_None_target_0"]

    on_torch = np.transpose(
        np.reshape(on_raw, [on_raw.shape[0], on_raw.shape[2], on_raw.shape[3], on_raw.shape[4]]),
        (3, 1, 2, 0),
    )[:, :, 0, :]
    off_torch = np.transpose(
        np.reshape(off_raw, [off_raw.shape[0], off_raw.shape[2], off_raw.shape[3], off_raw.shape[4]]),
        (3, 1, 2, 0),
    )[:, :, 0, :]
    noise_torch = np.transpose(noise_raw, (0, 3, 1, 2))[:, :, :, 0]

    on_single = on_torch[:, :, None, :sim_len].astype(np.float32)
    off_single = off_torch[:, :, None, :sim_len].astype(np.float32)
    noise_single = noise_torch[:, :, :sim_len, None].astype(np.float32)

    return {
        "on_torch": on_torch[:, :, :sim_len].astype(np.float32),
        "off_torch": off_torch[:, :, :sim_len].astype(np.float32),
        "noise_torch": noise_torch[:, :, :sim_len].astype(np.float32),
        "on_eprop_single": on_single,
        "off_eprop_single": off_single,
        "noise_eprop_single": noise_single,
        "rate_on": np.reshape(spks[on_key]["stimulus_0_rate"][:sim_len], (sim_len, batch_size)).astype(np.float32),
        "rate_off": np.reshape(spks[off_key]["stimulus_0_rate"][:sim_len], (sim_len, batch_size)).astype(np.float32),
        "rate_on_deriv": np.reshape(spks[on_key]["stimulus_0_rate_deriv"][:sim_len], (sim_len, batch_size)).astype(np.float32),
        "rate_off_deriv": np.reshape(spks[off_key]["stimulus_0_rate_deriv"][:sim_len], (sim_len, batch_size)).astype(np.float32),
    }


def bin_raster_batched(raster: torch.Tensor, bin_steps: int) -> torch.Tensor:
    batch, trials, time_steps = raster.shape
    n_bins = time_steps // bin_steps
    raster = raster[:, :, : n_bins * bin_steps]
    return raster.reshape(batch, trials, n_bins, bin_steps).sum(dim=(1, 3))


class BatchedSingleCellLIFBackprop(nn.Module):
    dt_ms = 0.1
    dt_s = 0.0001
    e_exc = 0.0
    e_inh = -80.0
    strf_latency_min_s = 0.001
    strf_latency_max_s = 0.05

    def __init__(self, params_2d: np.ndarray, surrogate_width: float = 0.5) -> None:
        super().__init__()
        params_2d = np.asarray(params_2d, dtype=np.float32)
        self.batch_size = params_2d.shape[1]
        self.surrogate_width = surrogate_width
        self.register_buffer("base_strf_gain", torch.as_tensor(np.maximum(params_2d[0], 1e-6)))
        self.register_buffer(
            "base_strf_latency",
            torch.as_tensor(np.clip(params_2d[1], self.strf_latency_min_s, self.strf_latency_max_s)),
        )
        self.raw_strf_gain = nn.Parameter(torch.as_tensor(inv_softplus_np(params_2d[0])))
        self.raw_strf_latency = nn.Parameter(
            torch.as_tensor(inv_sigmoid_range_np(params_2d[1], self.strf_latency_min_s, self.strf_latency_max_s))
        )
        self.raw_ron_g_inc = nn.Parameter(torch.as_tensor(inv_softplus_np(params_2d[2])))
        self.raw_gsyn = nn.ParameterDict(
            {
                "on_ron": nn.Parameter(torch.as_tensor(inv_softplus_np(params_2d[3]))),
                "off_ron": nn.Parameter(torch.as_tensor(inv_softplus_np(params_2d[4]))),
                "on_sonoff": nn.Parameter(torch.as_tensor(inv_softplus_np(params_2d[5]))),
                "off_sonoff": nn.Parameter(torch.as_tensor(inv_softplus_np(params_2d[6]))),
                "sonoff_ron": nn.Parameter(torch.as_tensor(inv_softplus_np(params_2d[7]))),
            }
        )

    def learned_params(self) -> np.ndarray:
        with torch.no_grad():
            out = np.zeros((8, self.batch_size), dtype=np.float32)
            out[0] = positive(self.raw_strf_gain).detach().cpu().numpy()
            out[1] = sigmoid_range(
                self.raw_strf_latency,
                self.strf_latency_min_s,
                self.strf_latency_max_s,
            ).detach().cpu().numpy()
            out[2] = positive(self.raw_ron_g_inc).detach().cpu().numpy()
            out[3] = positive(self.raw_gsyn["on_ron"]).detach().cpu().numpy()
            out[4] = positive(self.raw_gsyn["off_ron"]).detach().cpu().numpy()
            out[5] = positive(self.raw_gsyn["on_sonoff"]).detach().cpu().numpy()
            out[6] = positive(self.raw_gsyn["off_sonoff"]).detach().cpu().numpy()
            out[7] = positive(self.raw_gsyn["sonoff_ron"]).detach().cpu().numpy()
            return out

    def initial_state(self, batch: int, trials: int, device: torch.device, dtype: torch.dtype) -> Dict[str, object]:
        shape = (batch, trials)
        z = torch.zeros(shape, device=device, dtype=dtype)
        o = torch.ones(shape, device=device, dtype=dtype)
        state: Dict[str, object] = {
            "on_v": torch.full(shape, -65.0, device=device, dtype=dtype),
            "off_v": torch.full(shape, -65.0, device=device, dtype=dtype),
            "sonoff_v": torch.full(shape, -57.0, device=device, dtype=dtype),
            "ron_v": torch.full(shape, -65.0, device=device, dtype=dtype),
            "on_gad": z.clone(),
            "off_gad": z.clone(),
            "sonoff_gad": z.clone(),
            "ron_gad": z.clone(),
            "on_ref": z.clone(),
            "off_ref": z.clone(),
            "sonoff_ref": z.clone(),
            "ron_ref": z.clone(),
            "ron_noise_s": z.clone(),
            "ron_noise_x": z.clone(),
        }
        for name in ("on_ron", "off_ron", "on_sonoff", "off_sonoff", "sonoff_ron"):
            state[f"{name}_s"] = z.clone()
            state[f"{name}_x"] = z.clone()
            state[f"{name}_f"] = o.clone()
            state[f"{name}_p"] = o.clone()
            state[f"{name}_q"] = o.clone()
        for name, steps in {"on_ron": 10, "off_ron": 10, "on_sonoff": 30, "off_sonoff": 30, "sonoff_ron": 5}.items():
            state[f"{name}_delay"] = [z.clone() for _ in range(steps)]
        return state

    @staticmethod
    def detach_state(state: Dict[str, object]) -> Dict[str, object]:
        out: Dict[str, object] = {}
        for key, value in state.items():
            if isinstance(value, torch.Tensor):
                out[key] = value.detach()
            elif isinstance(value, list):
                out[key] = [item.detach() for item in value]
            else:
                out[key] = value
        return out

    @staticmethod
    def _fractional_shift_time(x: torch.Tensor, shift_steps: torch.Tensor) -> torch.Tensor:
        batch, trials, time_steps = x.shape
        time_idx = torch.arange(time_steps, device=x.device, dtype=x.dtype).view(1, time_steps)
        source = time_idx - shift_steps.view(batch, 1)
        idx0_float = torch.floor(source)
        frac = source - idx0_float
        idx0 = idx0_float.to(torch.long)
        idx1 = idx0 + 1
        valid0 = (idx0 >= 0) & (idx0 < time_steps)
        valid1 = (idx1 >= 0) & (idx1 < time_steps)
        idx0 = idx0.clamp(0, time_steps - 1)
        idx1 = idx1.clamp(0, time_steps - 1)
        idx0 = idx0[:, None, :].expand(batch, trials, time_steps)
        idx1 = idx1[:, None, :].expand(batch, trials, time_steps)
        y0 = torch.gather(x, 2, idx0) * valid0[:, None, :].to(x.dtype)
        y1 = torch.gather(x, 2, idx1) * valid1[:, None, :].to(x.dtype)
        frac = frac[:, None, :]
        return (1.0 - frac) * y0 + frac * y1

    def apply_strf_transform(self, on_input: torch.Tensor, off_input: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        batch = on_input.shape[0]
        strf_gain = positive(self.raw_strf_gain[:batch]).view(batch, 1, 1)
        strf_latency = sigmoid_range(
            self.raw_strf_latency[:batch],
            self.strf_latency_min_s,
            self.strf_latency_max_s,
        )
        base_gain = self.base_strf_gain[:batch].view(batch, 1, 1).clamp_min(1e-6)
        base_latency = self.base_strf_latency[:batch]
        gain_ratio = strf_gain / base_gain
        latency_delta_steps = (strf_latency - base_latency) / self.dt_s
        return (
            self._fractional_shift_time(on_input, latency_delta_steps) * gain_ratio,
            self._fractional_shift_time(off_input, latency_delta_steps) * gain_ratio,
        )

    def _spike_neuron(
        self,
        v: torch.Tensor,
        g_ad: torch.Tensor,
        ref: torch.Tensor,
        d_v: torch.Tensor,
        tau_ad: float,
        g_inc: torch.Tensor,
        v_thresh: float,
        v_reset: float,
        t_ref_ms: float,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        in_ref = ref > 0
        v_next = v + self.dt_ms * d_v
        g_next = g_ad + self.dt_ms * (-g_ad / tau_ad)
        spike = surrogate_spike(v_next - v_thresh, self.surrogate_width)
        spike = spike * (~in_ref).to(spike.dtype)
        spiked = spike.detach() > 0
        v_next = torch.where(spiked | in_ref, torch.as_tensor(v_reset, device=v.device, dtype=v.dtype), v_next)
        g_next = g_next + spike * g_inc
        ref_steps = int(round(t_ref_ms / self.dt_ms))
        ref_next = torch.clamp(ref - 1, min=0)
        ref_next = torch.where(spiked, torch.full_like(ref_next, ref_steps), ref_next)
        return v_next, g_next, ref_next, spike

    def _decay_synapse(self, state: Dict[str, object], name: str, tau_r: float, tau_d: float, tau_p: float, scale: float) -> None:
        s = state[f"{name}_s"]
        x = state[f"{name}_x"]
        f = state[f"{name}_f"]
        p = state[f"{name}_p"]
        state[f"{name}_s"] = s + self.dt_ms * ((scale * x - s) / tau_r)
        state[f"{name}_x"] = x + self.dt_ms * (-x / tau_d)
        state[f"{name}_f"] = f + self.dt_ms * ((1.0 - f) / 180.0)
        state[f"{name}_p"] = p + self.dt_ms * ((1.0 - p) / tau_p)

    def _delayed_event(self, state: Dict[str, object], name: str, spike: torch.Tensor) -> torch.Tensor:
        buffer = state[f"{name}_delay"]
        event = buffer.pop(0)
        buffer.append(spike)
        return event

    def _apply_synapse_event(self, state: Dict[str, object], name: str, event: torch.Tensor, f_p: float) -> None:
        x = state[f"{name}_x"]
        f = state[f"{name}_f"]
        p = state[f"{name}_p"]
        q = state[f"{name}_q"]
        state[f"{name}_x"] = x + q * event
        state[f"{name}_q"] = torch.where(event.detach() > 0, f * p, q)
        state[f"{name}_p"] = p * (1.0 - event * f_p)

    def forward(
        self,
        on_input: torch.Tensor,
        off_input: torch.Tensor,
        noise_input: torch.Tensor,
        bin_steps: int,
        state: Optional[Dict[str, object]] = None,
        return_trace: bool = False,
    ) -> Tuple[torch.Tensor, Dict[str, object], Optional[torch.Tensor]]:
        batch, trials, time_steps = on_input.shape
        if state is None:
            state = self.initial_state(batch, trials, on_input.device, on_input.dtype)
        on_input, off_input = self.apply_strf_transform(on_input, off_input)
        g = {name: positive(raw).view(batch, 1) for name, raw in self.raw_gsyn.items()}
        ron_g_inc = positive(self.raw_ron_g_inc).view(batch, 1)
        zero_g_inc = torch.zeros((batch, 1), device=on_input.device, dtype=on_input.dtype)
        input_g_inc = torch.full((batch, 1), 0.0003, device=on_input.device, dtype=on_input.dtype)

        psth_bins = []
        trace = [] if return_trace else None
        bin_accum = torch.zeros((batch, trials), device=on_input.device, dtype=on_input.dtype)
        for step in range(time_steps):
            on_t = on_input[:, :, step]
            off_t = off_input[:, :, step]
            noise_t = noise_input[:, :, step]
            self._decay_synapse(state, "on_ron", 0.7, 1.5, 30.0, 1.9481350796278847)
            self._decay_synapse(state, "off_ron", 0.7, 1.5, 30.0, 1.9481350796278847)
            self._decay_synapse(state, "on_sonoff", 0.1, 1.0, 80.0, 1.2915496650148839)
            self._decay_synapse(state, "off_sonoff", 0.1, 1.0, 80.0, 1.2915496650148839)
            self._decay_synapse(state, "sonoff_ron", 1.0, 4.5, 120.0, 1.5368523544529802)

            ron_noise_s = state["ron_noise_s"]
            ron_noise_x = state["ron_noise_x"]
            state["ron_noise_s"] = ron_noise_s + self.dt_ms * ((1.9481350796278847 * ron_noise_x - ron_noise_s) / 0.7)
            state["ron_noise_x"] = ron_noise_x + self.dt_ms * (-(ron_noise_x / 1.5) + noise_t / self.dt_ms)

            on_dv = ((-65.0 - state["on_v"]) - 200.0 * state["on_gad"] * (state["on_v"] + 80.0) - 200.0 * 0.17 * on_t * state["on_v"]) / 20.0
            off_dv = ((-65.0 - state["off_v"]) - 200.0 * state["off_gad"] * (state["off_v"] + 80.0) - 200.0 * 0.17 * off_t * state["off_v"]) / 20.0
            state["on_v"], state["on_gad"], state["on_ref"], on_spike = self._spike_neuron(state["on_v"], state["on_gad"], state["on_ref"], on_dv, 10.0, input_g_inc, -47.0, -54.0, 1.0)
            state["off_v"], state["off_gad"], state["off_ref"], off_spike = self._spike_neuron(state["off_v"], state["off_gad"], state["off_ref"], off_dv, 10.0, input_g_inc, -47.0, -54.0, 1.0)

            sonoff_syn = g["on_sonoff"] * state["on_sonoff_s"] * state["sonoff_v"] + g["off_sonoff"] * state["off_sonoff_s"] * state["sonoff_v"]
            sonoff_dv = ((-57.0 - state["sonoff_v"]) - 100.0 * state["sonoff_gad"] * (state["sonoff_v"] + 80.0) - 100.0 * sonoff_syn) / 10.0
            state["sonoff_v"], state["sonoff_gad"], state["sonoff_ref"], sonoff_spike = self._spike_neuron(state["sonoff_v"], state["sonoff_gad"], state["sonoff_ref"], sonoff_dv, 100.0, zero_g_inc, -47.0, -52.0, 0.5)

            ron_syn = (
                g["on_ron"] * state["on_ron_s"] * state["ron_v"]
                + g["off_ron"] * state["off_ron_s"] * state["ron_v"]
                + g["sonoff_ron"] * state["sonoff_ron_s"] * (state["ron_v"] + 80.0)
            )
            noise_current = -200.0 * 0.015 * ron_noise_s * state["ron_v"]
            ron_dv = ((-65.0 - state["ron_v"]) - 200.0 * state["ron_gad"] * (state["ron_v"] + 80.0) - 200.0 * ron_syn) / 20.0 + noise_current / 20.0
            state["ron_v"], state["ron_gad"], state["ron_ref"], ron_spike = self._spike_neuron(state["ron_v"], state["ron_gad"], state["ron_ref"], ron_dv, 100.0, ron_g_inc, -47.0, -54.0, 4.0)

            self._apply_synapse_event(state, "on_ron", self._delayed_event(state, "on_ron", on_spike), 0.1)
            self._apply_synapse_event(state, "off_ron", self._delayed_event(state, "off_ron", off_spike), 0.1)
            self._apply_synapse_event(state, "on_sonoff", self._delayed_event(state, "on_sonoff", on_spike), 0.2)
            self._apply_synapse_event(state, "off_sonoff", self._delayed_event(state, "off_sonoff", off_spike), 0.0)
            self._apply_synapse_event(state, "sonoff_ron", self._delayed_event(state, "sonoff_ron", sonoff_spike), 0.5)

            bin_accum = bin_accum + ron_spike
            if return_trace:
                trace.append(ron_spike)
            if (step + 1) % bin_steps == 0:
                psth_bins.append(bin_accum.sum(dim=1))
                bin_accum = torch.zeros_like(bin_accum)
        psth = torch.stack(psth_bins, dim=1) if psth_bins else torch.empty((batch, 0), device=on_input.device)
        trace_tensor = torch.stack(trace, dim=2) if return_trace and trace else None
        return psth, state, trace_tensor


def make_backprop_culprits(
    batch_size: int,
    trials: int,
    time_steps: int,
    bin_steps: int,
    known_tensor_bytes: int,
    param_and_optimizer_bytes: int,
    peak_allocated_mb: Optional[float],
) -> List[Dict[str, object]]:
    bins = time_steps // bin_steps
    binned_bytes = batch_size * bins * 4 * 3
    best_trace_bytes = batch_size * trials * time_steps * 4
    peak_bytes = (peak_allocated_mb or 0.0) * 1024**2
    graph_est = max(0.0, peak_bytes - known_tensor_bytes - param_and_optimizer_bytes - binned_bytes)
    culprits = [
        ("autograd saved timestep graph and intermediates", graph_est, "measured residual from PyTorch peak allocation"),
        ("input/target/noise tensors", known_tensor_bytes, "measured tensor sizes"),
        ("optimizer state, parameters, and gradients", param_and_optimizer_bytes, "measured/estimated tensor sizes"),
        ("binned PSTH/loss tensors", binned_bytes, "estimated from bin count"),
        ("best-output trace snapshot", best_trace_bytes, "estimated no_grad trace copy"),
    ]
    return [{"name": n, "mb": mb(v), "source": s} for n, v, s in culprits if v > 0]


def make_eprop_culprits(batch_size: int, trials: int, time_steps: int, backend: str) -> List[Dict[str, object]]:
    array_backend = "CuPy/GPU" if backend == "cuda" else "NumPy/CPU"
    num_cells = 1
    channels = 1
    output_holders = 4 * num_cells * batch_size * trials * channels * time_steps
    gpu_inputs = 4 * (3 * num_cells * batch_size * trials * channels * time_steps + num_cells * trials * time_steps * channels)
    eligibility = 8 * (5 * num_cells * batch_size * trials * channels * 100 + 8 * num_cells * batch_size * trials * channels)
    state = 8 * (
        4 * num_cells * batch_size * trials * channels * (2 + 2 + 5 + 1)
        + 5 * num_cells * batch_size * trials * channels * (2 * 5)
    )
    rate_arrays = 4 * 4 * time_steps * num_cells * batch_size
    culprits = [
        ("output/PSC/loss spike holders", output_holders, f"estimated generated {array_backend} int8 holders"),
        ("input, noise, and target arrays", gpu_inputs, f"estimated {array_backend} float32 inputs/data"),
        ("eligibility and gradient accumulators", eligibility, f"estimated generated {array_backend} float64 accumulators"),
        ("neuron and synapse rolling state", state, f"estimated generated {array_backend} float64 state"),
        ("STRF rate and derivative arrays", rate_arrays, "estimated float32 rate arrays"),
    ]
    return [{"name": n, "mb": mb(v), "source": s} for n, v, s in culprits if v > 0]


def capture_backprop_trace(
    model: BatchedSingleCellLIFBackprop,
    on_t: torch.Tensor,
    off_t: torch.Tensor,
    noise_t: torch.Tensor,
    bin_steps: int,
) -> np.ndarray:
    was_training = model.training
    model.eval()
    with torch.no_grad():
        _, _, trace = model(on_t, off_t, noise_t, bin_steps, return_trace=True)
    if was_training:
        model.train()
    return trace.detach().cpu().numpy().astype(np.float32)


def run_backprop_method(
    method: str,
    params_2d: np.ndarray,
    inputs: Dict[str, np.ndarray],
    target_raster: np.ndarray,
    epochs: int,
    bin_steps: int,
    device: torch.device,
    lr: float,
    tbptt_steps: Optional[int],
) -> Tuple[BenchResult, Dict[str, np.ndarray]]:
    batch_size = params_2d.shape[1]
    on_t = torch.as_tensor(inputs["on_torch"], dtype=torch.float32, device=device)
    off_t = torch.as_tensor(inputs["off_torch"], dtype=torch.float32, device=device)
    noise_t = torch.as_tensor(inputs["noise_torch"], dtype=torch.float32, device=device)
    target = np.repeat(target_raster[None, :, : on_t.shape[-1]], batch_size, axis=0)
    target_t = torch.as_tensor(target, dtype=torch.float32, device=device)
    target_bins = bin_raster_batched(target_t, bin_steps)

    model = BatchedSingleCellLIFBackprop(params_2d).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    known_bytes = tensor_bytes(on_t, off_t, noise_t, target_t)
    param_bytes = sum(p.numel() * p.element_size() for p in model.parameters()) * 4
    phase = PhaseProfile()
    best_loss = np.inf
    best_epoch = 0
    best_batch_id = 1
    best_output = capture_backprop_trace(model, on_t, off_t, noise_t, bin_steps)
    initial_output = best_output.copy()
    loss_history = []

    rss_start = current_rss_mb()
    reset_torch_memory(device)
    start_all = time.perf_counter()
    for epoch in range(epochs):
        state = None
        epoch_loss_parts = []
        if tbptt_steps is None:
            z0 = time.perf_counter()
            optimizer.zero_grad(set_to_none=True)
            sync_torch(device)
            phase.zero_grad_s += time.perf_counter() - z0
            f0 = time.perf_counter()
            pred_bins, _, _ = model(on_t, off_t, noise_t, bin_steps)
            loss_per_batch = ((pred_bins - target_bins) ** 2).sum(dim=1)
            loss = loss_per_batch.mean()
            sync_torch(device)
            phase.forward_s += time.perf_counter() - f0
            b0 = time.perf_counter()
            loss.backward()
            sync_torch(device)
            phase.backward_s += time.perf_counter() - b0
            s0 = time.perf_counter()
            optimizer.step()
            sync_torch(device)
            phase.optimizer_s += time.perf_counter() - s0
            epoch_loss_np = loss_per_batch.detach().cpu().numpy()
        else:
            usable = (on_t.shape[-1] // tbptt_steps) * tbptt_steps
            for start_idx in range(0, usable, tbptt_steps):
                stop_idx = start_idx + tbptt_steps
                z0 = time.perf_counter()
                optimizer.zero_grad(set_to_none=True)
                sync_torch(device)
                phase.zero_grad_s += time.perf_counter() - z0
                f0 = time.perf_counter()
                pred_bins, state, _ = model(
                    on_t[:, :, start_idx:stop_idx],
                    off_t[:, :, start_idx:stop_idx],
                    noise_t[:, :, start_idx:stop_idx],
                    bin_steps,
                    state=state,
                )
                target_chunk_bins = bin_raster_batched(target_t[:, :, start_idx:stop_idx], bin_steps)
                loss_per_batch = ((pred_bins - target_chunk_bins) ** 2).sum(dim=1)
                loss = loss_per_batch.mean()
                sync_torch(device)
                phase.forward_s += time.perf_counter() - f0
                b0 = time.perf_counter()
                loss.backward()
                sync_torch(device)
                phase.backward_s += time.perf_counter() - b0
                s0 = time.perf_counter()
                optimizer.step()
                sync_torch(device)
                phase.optimizer_s += time.perf_counter() - s0
                epoch_loss_parts.append(loss_per_batch.detach().cpu().numpy())
                state = model.detach_state(state)
            epoch_loss_np = np.mean(np.stack(epoch_loss_parts, axis=0), axis=0)

        loss_history.append(epoch_loss_np)
        idx = int(np.argmin(epoch_loss_np))
        if float(epoch_loss_np[idx]) < best_loss:
            t0 = time.perf_counter()
            trace = capture_backprop_trace(model, on_t, off_t, noise_t, bin_steps)
            sync_torch(device)
            phase.best_trace_s += time.perf_counter() - t0
            best_output = trace
            best_loss = float(epoch_loss_np[idx])
            best_batch_id = idx + 1
            best_epoch = epoch + 1

    phase.total_s = time.perf_counter() - start_all
    rss_end = current_rss_mb()
    peak_allocated, peak_reserved = torch_memory_mb(device)
    culprits = make_backprop_culprits(
        batch_size=batch_size,
        trials=on_t.shape[1],
        time_steps=on_t.shape[2],
        bin_steps=bin_steps,
        known_tensor_bytes=known_bytes,
        param_and_optimizer_bytes=param_bytes,
        peak_allocated_mb=peak_allocated,
    )
    phase_memory = {
        "forward_graph_and_intermediates_mb": culprits[0]["mb"] if culprits else 0.0,
        "backward_extra_peak_mb": 0.0,
        "inputs_targets_mb": mb(known_bytes),
        "optimizer_params_grads_mb": mb(param_bytes),
    }
    result = BenchResult(
        test="",
        method=method,
        batch_size=batch_size,
        epochs=epochs,
        final_loss=float(np.mean(loss_history[-1])),
        best_loss=best_loss,
        best_batch_id=best_batch_id,
        best_epoch=best_epoch,
        wall_time_s=phase.total_s,
        peak_cuda_allocated_mb=peak_allocated,
        peak_cuda_reserved_mb=peak_reserved,
        cupy_pool_used_mb=None,
        cupy_pool_reserved_mb=None,
        rss_start_mb=rss_start,
        rss_end_mb=rss_end,
        rss_delta_mb=None if rss_start is None or rss_end is None else rss_end - rss_start,
        phase_profile=asdict(phase),
        memory_culprits=culprits,
        phase_memory=phase_memory,
    )
    artifacts = {
        "initial_output": initial_output,
        "best_output": best_output,
        "loss_history": np.asarray(loss_history),
        "learned_params": model.learned_params(),
    }
    return result, artifacts


def patch_single_cell_eprop_body(body: str) -> str:
    """Fix single-cell generated binned-loss reductions for benchmark shapes."""
    body = body.replace(
        "            sim_bin = np.sum(np.squeeze(np.sum(ROn_spikes_holder[:,:,:,timestep-loss_bin_width:timestep], axis=-1)),axis=-1)",
        "            sim_bin = np.sum(ROn_spikes_holder[:,:,:,timestep-loss_bin_width:timestep], axis=(1,2,3))",
    )
    body = body.replace(
        "            data_bin = np.sum(np.sum(np.squeeze(data[:,timestep-loss_bin_width:timestep,:], axis=-1)),axis=-1)",
        "            data_bin = np.sum(data[:,timestep-loss_bin_width:timestep,:])",
    )
    return body


def make_single_cell_cupy_body(body: str) -> str:
    body = patch_single_cell_eprop_body(body)
    body = body.replace("np.", "cp.")
    transfer_marker = "\n\n    #Declare Variables\n"
    transfers = (
        "\n\n    #Transfer inputs to GPU\n"
        "    on_input = cp.asarray(on_input)\n"
        "    off_input = cp.asarray(off_input)\n"
        "    noise_token = cp.asarray(noise_token)\n"
        "    rate_on = cp.asarray(rate_on)\n"
        "    rate_off = cp.asarray(rate_off)\n"
        "    rate_on_deriv = cp.asarray(rate_on_deriv)\n"
        "    rate_off_deriv = cp.asarray(rate_off_deriv)\n"
        "    data = cp.asarray(data)\n"
        "    p = cp.asarray(p)\n"
    )
    if transfer_marker not in body:
        raise RuntimeError("Could not find generated single-cell variable declaration marker.")
    body = body.replace(transfer_marker, transfers + transfer_marker, 1)
    return_line = "    return ROn_spikes_holder, grads, On_SOnOff_PSC_s_holder, Off_SOnOff_PSC_s_holder, losses_holder"
    cupy_return = (
        "    cp.cuda.Stream.null.synchronize()\n"
        "    return cp.asnumpy(ROn_spikes_holder), cp.asnumpy(grads), "
        "cp.asnumpy(On_SOnOff_PSC_s_holder), cp.asnumpy(Off_SOnOff_PSC_s_holder), "
        "cp.asnumpy(losses_holder)"
    )
    if return_line not in body:
        raise RuntimeError("Could not find generated single-cell return line.")
    return body.replace(return_line, cupy_return, 1)


def compile_single_cell_eprop_solver(batch_size: int, sim_len: int, results_dir: Path, backend: str):
    with eprop_import_path():
        import declarations
        import set_options
        from BuildFile import Forwards_Method

        opts = set_options.options()
        opts["N_batch"] = batch_size
        opts["sim_len"] = sim_len
        arch = declarations.Declare_Architecture(opts)
        body = Forwards_Method.Euler_Compiler(arch[0], arch[1], arch[2], opts)

    if backend == "cuda":
        body = make_single_cell_cupy_body(body)
        imports = "import numpy as np\nimport cupy as cp\nfrom BuildFile import calculate_loss_eprop\n"
    elif backend == "cpu":
        body = patch_single_cell_eprop_body(body)
        imports = "import numpy as np\nfrom BuildFile import calculate_loss_eprop\n"
    else:
        raise ValueError(f"Unknown E-prop backend: {backend}")

    solver_path = results_dir / f"generated_single_cell_eprop_{backend}_bs{batch_size}_t{sim_len}.py"
    solver_text = (
        imports
        + "def solve_run(on_input,off_input,noise_token,rate_on,rate_off,rate_on_deriv,rate_off_deriv,data,p):\n"
        + body
    )
    solver_path.write_text(solver_text, encoding="utf-8")
    module_name = f"generated_single_cell_eprop_{backend}_bs{batch_size}_t{sim_len}_{int(time.time() * 1000)}"
    with eprop_import_path():
        spec = importlib.util.spec_from_file_location(module_name, solver_path)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Could not import generated solver {solver_path}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return module


def eprop_adam_update(p: np.ndarray, m: np.ndarray, v: np.ndarray, t: int, grads: np.ndarray, lr: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    beta1, beta2, eps = 0.5, 0.9995, 1e-6
    t += 1
    m = beta1 * m + (1 - beta1) * grads
    v = beta2 * v + (1 - beta2) * (grads**2)
    m_hat = m / (1 - beta1**t)
    v_hat = v / (1 - beta2**t)
    p = p - lr * m_hat / (np.sqrt(v_hat) + eps)
    p[0:2, :] = np.where(p[0:2, :] < 0.001, 0.001, p[0:2, :])
    p[1, :] = np.where(p[1, :] > 0.05, 0.05, p[1, :])
    p[2:8, :] = np.where(p[2:8, :] < 0, 0, p[2:8, :])
    return p, m, v, t


def run_eprop_method(
    method: str,
    backend: str,
    params_2d: np.ndarray,
    inputs: Dict[str, np.ndarray],
    target_raster: np.ndarray,
    epochs: int,
    bin_steps: int,
    lr: float,
    results_dir: Path,
) -> Tuple[BenchResult, Dict[str, np.ndarray]]:
    batch_size = params_2d.shape[1]
    sim_len = inputs["on_eprop_single"].shape[-1]
    solver = compile_single_cell_eprop_solver(batch_size, sim_len, results_dir, backend)
    with eprop_import_path():
        from BuildFile import calculate_loss

    data = target_raster[:, :sim_len, None].astype(np.float32)
    p = params_2d.astype(np.float32).copy()
    m = np.zeros_like(p)
    v = np.zeros_like(p)
    t = 0
    loss_history = []
    best_loss = np.inf
    best_batch_id = 1
    best_epoch = 0
    best_output = None
    initial_output = None
    phase = PhaseProfile()

    if backend == "cuda":
        reset_cupy_memory()
    rss_start = current_rss_mb()
    start_all = time.perf_counter()
    for epoch in range(epochs):
        f0 = time.perf_counter()
        output, grads, _, _, _ = solver.solve_run(
            inputs["on_eprop_single"],
            inputs["off_eprop_single"],
            inputs["noise_eprop_single"],
            inputs["rate_on"],
            inputs["rate_off"],
            inputs["rate_on_deriv"],
            inputs["rate_off_deriv"],
            data,
            p,
        )
        if backend == "cuda":
            cp.cuda.Stream.null.synchronize()
        phase.forward_s += time.perf_counter() - f0
        if initial_output is None:
            initial_output = output.copy()
        _, loss, _ = calculate_loss.calculate(output, grads, data)
        loss_per_batch = np.asarray(loss[1], dtype=np.float32).reshape(batch_size)
        loss_history.append(loss_per_batch)
        idx = int(np.argmin(loss_per_batch))
        if float(loss_per_batch[idx]) < best_loss:
            best_loss = float(loss_per_batch[idx])
            best_batch_id = idx + 1
            best_epoch = epoch + 1
            best_output = output[idx, :, 0, :].astype(np.float32)
        b0 = time.perf_counter()
        grads2d = np.asarray(grads, dtype=np.float32).reshape(8, batch_size)
        p, m, v, t = eprop_adam_update(p, m, v, t, grads2d, lr)
        phase.optimizer_s += time.perf_counter() - b0

    phase.total_s = time.perf_counter() - start_all
    rss_end = current_rss_mb()
    if backend == "cuda":
        used_mb, reserved_mb = cupy_pool_mb()
    else:
        used_mb, reserved_mb = None, None
    if best_output is None:
        best_output = np.zeros((target_raster.shape[0], sim_len), dtype=np.float32)
    if initial_output is None:
        initial_output = np.zeros((batch_size, target_raster.shape[0], 1, sim_len), dtype=np.float32)
    culprits = make_eprop_culprits(batch_size, target_raster.shape[0], sim_len, backend)
    phase_memory = {
        "forward_state_and_output_holders_mb": culprits[0]["mb"],
        "eligibility_gradient_holders_mb": culprits[2]["mb"],
        "inputs_targets_mb": culprits[1]["mb"],
    }
    result = BenchResult(
        test="",
        method=method,
        batch_size=batch_size,
        epochs=epochs,
        final_loss=float(np.mean(loss_history[-1])),
        best_loss=best_loss,
        best_batch_id=best_batch_id,
        best_epoch=best_epoch,
        wall_time_s=phase.total_s,
        peak_cuda_allocated_mb=None,
        peak_cuda_reserved_mb=None,
        cupy_pool_used_mb=used_mb,
        cupy_pool_reserved_mb=reserved_mb,
        rss_start_mb=rss_start,
        rss_end_mb=rss_end,
        rss_delta_mb=None if rss_start is None or rss_end is None else rss_end - rss_start,
        phase_profile=asdict(phase),
        memory_culprits=culprits,
        phase_memory=phase_memory,
    )
    artifacts = {
        "initial_output": initial_output[best_batch_id - 1, :, 0, :].astype(np.float32),
        "best_output": best_output,
        "loss_history": np.asarray(loss_history),
        "learned_params": p,
    }
    return result, artifacts


def write_csv(path: Path, rows: List[Dict[str, object]], fields: List[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def psth(raster: np.ndarray, bin_steps: int) -> np.ndarray:
    if raster.ndim == 3:
        raster = raster[0]
    n_bins = raster.shape[-1] // bin_steps
    raster = raster[:, : n_bins * bin_steps]
    return raster.reshape(raster.shape[0], n_bins, bin_steps).sum(axis=(0, 2))


def plot_speed_outputs(path: Path, target: np.ndarray, artifacts_by_method: Dict[str, Dict[str, np.ndarray]], bin_steps: int) -> None:
    if plt is None:
        return
    methods = list(artifacts_by_method.keys())
    fig, axes = plt.subplots(len(methods), 2, figsize=(10, 3.0 * len(methods)), squeeze=False)
    target_psth = psth(target, bin_steps)
    for row, method in enumerate(methods):
        art = artifacts_by_method[method]
        axes[row, 0].plot(target_psth, color="black", lw=1.5, label="data")
        axes[row, 0].plot(psth(art["initial_output"], bin_steps), color="#8a8f98", lw=1.2, label="start")
        axes[row, 0].set_title(f"{method}: starting output")
        axes[row, 1].plot(target_psth, color="black", lw=1.5, label="data")
        axes[row, 1].plot(psth(art["best_output"], bin_steps), color="#2563eb", lw=1.2, label="best")
        axes[row, 1].set_title(f"{method}: best after optimization")
        for ax in axes[row]:
            ax.set_xlabel("PSTH bin")
            ax.set_ylabel("spikes / bin")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_report(results_dir: Path, results: List[BenchResult], args: argparse.Namespace) -> None:
    lines = [
        "# E-prop vs Backprop Benchmark Report",
        "",
        f"Cell: {args.cell_id}",
        f"LongRunResults source: `{args.long_run_mat}`",
        f"Sim length: {args.sim_len} steps",
        f"Bin steps: {args.bin_steps}",
        "",
        "## Memory Summary",
        "",
        "| Method | Batch | Epochs | Wall s | RSS delta MB | Torch peak allocated MB | Torch peak reserved MB | CuPy pool reserved MB | Best loss |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in results:
        if r.test != "memory":
            continue
        lines.append(
            f"| {r.method} | {r.batch_size} | {r.epochs} | {r.wall_time_s:.2f} | "
            f"{'' if r.rss_delta_mb is None else f'{r.rss_delta_mb:.1f}'} | "
            f"{'' if r.peak_cuda_allocated_mb is None else f'{r.peak_cuda_allocated_mb:.1f}'} | "
            f"{'' if r.peak_cuda_reserved_mb is None else f'{r.peak_cuda_reserved_mb:.1f}'} | "
            f"{'' if r.cupy_pool_reserved_mb is None else f'{r.cupy_pool_reserved_mb:.1f}'} | "
            f"{r.best_loss:.4f} |"
        )
    lines += [
        "",
        "## Speed Summary",
        "",
        "| Method | Batch | Epochs | Wall s | s / epoch | Best loss | Best batch | Best epoch |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in results:
        if r.test != "speed":
            continue
        lines.append(
            f"| {r.method} | {r.batch_size} | {r.epochs} | {r.wall_time_s:.2f} | "
            f"{r.wall_time_s / max(r.epochs, 1):.2f} | {r.best_loss:.4f} | "
            f"{r.best_batch_id} | {r.best_epoch} |"
        )
    lines += [
        "",
        "## Memory Culprits",
        "",
        "These culprits are broad, non-overlapping categories. PyTorch peak memory is measured through "
        "`torch.cuda.max_memory_allocated`; CuPy/E-prop CUDA memory is pool-measured plus generated-array estimates. "
        "CPU memory uses process RSS delta when `psutil` is available.",
        "",
    ]
    for r in results:
        lines.append(f"### {r.test}: {r.method}, batch {r.batch_size}")
        for item in r.memory_culprits:
            lines.append(f"- {item['name']}: {item['mb']:.2f} MB ({item['source']})")
        lines.append("")
    lines += [
        "## Interpretation Checks",
        "",
        "- Full backprop stores a timestep-by-timestep autograd graph until `loss.backward()` finishes.",
        "- Truncated backprop should lower peak graph memory roughly in proportion to the chunk length, "
        "but it changes the gradient because state is detached between chunks.",
        "- E-prop does not store a full autograd graph. Its persistent memory should be dominated by "
        "generated holders, eligibility accumulators, and input/output arrays.",
        "- `eprop_cpu` is generated from `E-prop/BuildFile/Forwards_Method.py`, the same single-cell path used by "
        "`run_main_integrated.py`. `eprop_cupy` is generated by converting that same single-cell body to CuPy.",
        "- The speed plot `speed_fit_outputs.png` compares data PSTH with starting and best outputs.",
    ]
    (results_dir / "benchmark_report.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-id", type=int, default=7)
    parser.add_argument("--long-run-mat", type=Path, default=make_default_long_run_path())
    parser.add_argument("--oliver-data-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=Path(__file__).resolve().parent / "benchmark_results")
    parser.add_argument("--batch-sizes", type=int, nargs="+", default=[5, 20])
    parser.add_argument("--speed-batch-size", type=int, default=5)
    parser.add_argument("--speed-epochs", type=int, default=10)
    parser.add_argument("--sim-len", type=int, default=29801)
    parser.add_argument("--dt-ms", type=float, default=0.1)
    parser.add_argument("--bin-steps", type=int, default=200)
    parser.add_argument("--tbptt-steps", type=int, default=2000)
    parser.add_argument("--lr", type=float, default=5e-1)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--backprop-devices", choices=["cpu", "cuda"], nargs="+", default=["cpu", "cuda"])
    parser.add_argument("--skip-eprop-cpu", action="store_true")
    parser.add_argument("--skip-eprop-cupy", action="store_true")
    parser.add_argument("--skip-memory", action="store_true")
    parser.add_argument("--skip-speed", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    backprop_devices: List[torch.device] = []
    for name in args.backprop_devices:
        if name == "cuda" and not torch.cuda.is_available():
            print("CUDA requested for backprop, but torch.cuda.is_available() is false. Skipping backprop_cuda.")
            continue
        backprop_devices.append(torch.device(name))
    if not backprop_devices:
        raise SystemExit("No usable backprop devices selected.")

    eprop_backends: List[Tuple[str, str]] = []
    if not args.skip_eprop_cpu:
        eprop_backends.append(("eprop_cpu", "cpu"))
    if not args.skip_eprop_cupy:
        if cp is None:
            raise SystemExit("CuPy is not installed, so eprop_cupy cannot run. Use --skip-eprop-cupy to run CPU-only E-prop.")
        eprop_backends.append(("eprop_cupy", "cuda"))

    from run_single_cell_backprop import default_oliver_data_dir

    oliver_dir = args.oliver_data_dir or default_oliver_data_dir()
    target_raster, spontaneous_fr, focus = load_cell_raster(args.cell_id, oliver_dir, args.dt_ms, args.sim_len, "peak")
    target_raster = target_raster[:, : args.sim_len].astype(np.float32)
    print(f"Benchmarking cell {args.cell_id}, focus column {focus + 1}, FR proxy {spontaneous_fr:.3f} Hz")
    print(f"Backprop devices: {', '.join(str(d) for d in backprop_devices)}")
    print(f"E-prop backends: {', '.join(name for name, _ in eprop_backends)}")

    all_results: List[BenchResult] = []
    speed_artifacts: Dict[str, Dict[str, np.ndarray]] = {}

    if not args.skip_memory:
        for batch_size in args.batch_sizes:
            params = load_cell_params(args.long_run_mat, args.cell_id, batch_size, perturb=False)
            inputs = generate_inputs(params, batch_size, spontaneous_fr, args.sim_len, args.seed)
            for device in backprop_devices:
                for base_method, tbptt in [("backprop", None), ("truncated_backprop", args.tbptt_steps)]:
                    method = f"{base_method}_{device.type}"
                    print(f"Memory test: {method}, batch {batch_size}")
                    result, _ = run_backprop_method(method, params, inputs, target_raster, 1, args.bin_steps, device, args.lr, tbptt)
                    result.test = "memory"
                    all_results.append(result)
            for method, backend in eprop_backends:
                print(f"Memory test: {method}, batch {batch_size}")
                result, _ = run_eprop_method(method, backend, params, inputs, target_raster, 1, args.bin_steps, args.lr, args.results_dir)
                result.test = "memory"
                all_results.append(result)

    if not args.skip_speed:
        batch_size = args.speed_batch_size
        params = load_cell_params(args.long_run_mat, args.cell_id, batch_size, perturb=True)
        inputs = generate_inputs(params, batch_size, spontaneous_fr, args.sim_len, args.seed + 100)
        for device in backprop_devices:
            for base_method, tbptt in [("backprop", None), ("truncated_backprop", args.tbptt_steps)]:
                method = f"{base_method}_{device.type}"
                print(f"Speed test: {method}, batch {batch_size}, epochs {args.speed_epochs}")
                result, artifacts = run_backprop_method(method, params, inputs, target_raster, args.speed_epochs, args.bin_steps, device, args.lr, tbptt)
                result.test = "speed"
                all_results.append(result)
                speed_artifacts[method] = artifacts
        for method, backend in eprop_backends:
            print(f"Speed test: {method}, batch {batch_size}, epochs {args.speed_epochs}")
            result, artifacts = run_eprop_method(method, backend, params, inputs, target_raster, args.speed_epochs, args.bin_steps, args.lr, args.results_dir)
            result.test = "speed"
            all_results.append(result)
            speed_artifacts[method] = artifacts

    json_payload = [asdict(r) for r in all_results]
    (args.results_dir / "benchmark_results.json").write_text(json.dumps(json_payload, indent=2), encoding="utf-8")
    summary_rows = [
        {
            "test": r.test,
            "method": r.method,
            "batch_size": r.batch_size,
            "epochs": r.epochs,
            "wall_time_s": r.wall_time_s,
            "seconds_per_epoch": r.wall_time_s / max(r.epochs, 1),
            "best_loss": r.best_loss,
            "best_batch_id": r.best_batch_id,
            "best_epoch": r.best_epoch,
            "torch_peak_allocated_mb": r.peak_cuda_allocated_mb,
            "torch_peak_reserved_mb": r.peak_cuda_reserved_mb,
            "cupy_pool_reserved_mb": r.cupy_pool_reserved_mb,
            "rss_delta_mb": r.rss_delta_mb,
        }
        for r in all_results
    ]
    fields = list(summary_rows[0].keys()) if summary_rows else []
    if fields:
        write_csv(args.results_dir / "benchmark_summary.csv", summary_rows, fields)
        write_csv(args.results_dir / "memory_summary.csv", [r for r in summary_rows if r["test"] == "memory"], fields)
        write_csv(args.results_dir / "speed_summary.csv", [r for r in summary_rows if r["test"] == "speed"], fields)
    if speed_artifacts:
        plot_speed_outputs(args.results_dir / "speed_fit_outputs.png", target_raster, speed_artifacts, args.bin_steps)
    write_report(args.results_dir, all_results, args)
    print(f"Benchmark complete. See {args.results_dir}")


if __name__ == "__main__":
    main()
