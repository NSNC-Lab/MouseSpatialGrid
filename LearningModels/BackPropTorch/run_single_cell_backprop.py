"""
Single-cell PyTorch backprop prototype for the MouseSpatialGrid E-prop model.

This is intentionally side-by-side with the existing E-prop code. It reuses the
current STRF/input generation path, then trains a single-cell LIF network with
PyTorch autograd and surrogate-gradient spikes.
"""

from __future__ import annotations

import argparse
import contextlib
import json
import math
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
except ModuleNotFoundError as exc:
    raise SystemExit(
        "PyTorch is required for this runner. Activate an environment with torch "
        "installed, then rerun this script."
    ) from exc

try:
    from scipy.io import loadmat, savemat
except ModuleNotFoundError as exc:
    raise SystemExit(
        "scipy is required to read the Oliver .mat file and save results. "
        "Install scipy in the same environment as torch."
    ) from exc


SCRIPT_DIR = Path(__file__).resolve().parent
LEARNING_MODELS_DIR = SCRIPT_DIR.parent
EPROP_DIR = LEARNING_MODELS_DIR / "E-prop"


@dataclass
class InitialParams:
    strf_gain: float
    strf_latency: float
    ron_g_inc: float
    on_ron_gsyn: float
    off_ron_gsyn: float
    on_sonoff_gsyn: float
    off_sonoff_gsyn: float
    sonoff_ron_gsyn: float

    def as_eprop_array(self, batch_size: int) -> np.ndarray:
        p = np.zeros((8, batch_size), dtype=np.float32)
        p[0, :] = self.strf_gain
        p[1, :] = self.strf_latency
        p[2, :] = self.ron_g_inc
        p[3, :] = self.on_ron_gsyn
        p[4, :] = self.off_ron_gsyn
        p[5, :] = self.on_sonoff_gsyn
        p[6, :] = self.off_sonoff_gsyn
        p[7, :] = self.sonoff_ron_gsyn
        return p


def make_initial_params(seed: int) -> InitialParams:
    rng = np.random.default_rng(seed)
    return InitialParams(
        strf_gain=float(rng.uniform(0.1, 0.2)),
        strf_latency=float(rng.uniform(0.02, 0.03)),
        ron_g_inc=float(rng.uniform(0.0, 1e-5)),
        on_ron_gsyn=float(rng.uniform(0.02, 0.2)),
        off_ron_gsyn=float(rng.uniform(0.02, 0.2)),
        on_sonoff_gsyn=float(rng.uniform(0.0, 1e-4)),
        off_sonoff_gsyn=float(rng.uniform(0.0, 1e-4)),
        sonoff_ron_gsyn=float(rng.uniform(0.0, 1e-4)),
    )


def inv_softplus(value: float) -> float:
    value = max(float(value), 1e-10)
    return math.log(math.expm1(value))


def positive(raw: torch.Tensor) -> torch.Tensor:
    return F.softplus(raw) + 1e-10


def surrogate_spike(v_minus_threshold: torch.Tensor, width: float) -> torch.Tensor:
    hard = (v_minus_threshold >= 0).to(v_minus_threshold.dtype)
    soft = torch.sigmoid(v_minus_threshold / width)
    return hard.detach() - soft.detach() + soft


@contextlib.contextmanager
def clean_argv() -> Iterable[None]:
    old_argv = sys.argv[:]
    sys.argv = [old_argv[0]]
    try:
        yield
    finally:
        sys.argv = old_argv


def get_focus_index(tuning_value: object) -> int:
    text = str(tuning_value).lower()
    if "contra" in text:
        return 0
    if "45" in text:
        return 1
    if "center" in text:
        return 2
    if "ipsi" in text:
        return 3
    return 0


def default_oliver_data_dir() -> Path:
    env_dir = os.environ.get("OLIVER_DATA_DIR")
    if env_dir:
        return Path(env_dir)

    candidates = [
        Path(r"C:\Users\ipboy\Documents\GitHub\ModelingEffort\Multi-Channel\Plotting\OliverDataPlotting"),
        Path(r"C:\Users\ipboy\Documents\Github\ModelingEffort\Multi-Channel\Plotting\OliverDataPlotting"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def load_cell_raster(
    cell_id: int,
    data_dir: Path,
    dt_ms: float,
    sim_len: int,
    focus_mode: str,
) -> Tuple[np.ndarray, float, int]:
    mat_path = data_dir / "all_units_info_with_polished_criteria_modified_perf.mat"
    if not mat_path.exists():
        raise FileNotFoundError(
            f"Could not find {mat_path}. Set OLIVER_DATA_DIR or pass --oliver-data-dir."
        )

    mat = loadmat(
        mat_path,
        variable_names=["all_data"],
        squeeze_me=True,
        struct_as_record=False,
    )
    all_data = mat["all_data"]
    unit = all_data[cell_id - 1]
    spike_times = unit.ctrl_tar1_timestamps

    focus = 0
    if focus_mode == "peak":
        focus = get_focus_index(getattr(unit, "tuning_type", ""))

    if isinstance(spike_times, np.ndarray) and spike_times.ndim == 2:
        focus = min(focus, spike_times.shape[1] - 1)
        spike_times = spike_times[:, focus]

    raster = np.zeros((10, sim_len), dtype=np.float32)
    pre_spikes = []
    stim_end_s = (sim_len - 1) * dt_ms / 1000.0
    samples_per_second = 1000.0 / dt_ms

    for trial in range(10):
        times = np.asarray(spike_times[trial]).squeeze()
        times = np.atleast_1d(times).astype(np.float64)
        times = times[np.isfinite(times)]
        pre_spikes.append(times[times < 0])

        stim_mask = (times >= 0) & (times < stim_end_s)
        indices = np.round(times[stim_mask] * samples_per_second).astype(np.int64)
        indices = indices[(indices >= 0) & (indices < sim_len)]
        raster[trial, indices] = 1.0

    pre_concat = np.concatenate(pre_spikes) if pre_spikes else np.array([])
    spontaneous_fr = float(pre_concat.size / 10.0)
    return raster, spontaneous_fr, focus


def trim_or_pad(x: np.ndarray, target_len: int, axis: int = -1) -> np.ndarray:
    axis = axis % x.ndim
    current = x.shape[axis]
    if current == target_len:
        return x
    if current > target_len:
        sl = [slice(None)] * x.ndim
        sl[axis] = slice(0, target_len)
        return x[tuple(sl)]

    pad_width = [(0, 0)] * x.ndim
    pad_width[axis] = (0, target_len - current)
    return np.pad(x, pad_width, mode="constant")


def generate_eprop_inputs(
    spontaneous_fr: float,
    params: InitialParams,
    sim_len: int,
    input_batch_size: int,
    seed: int,
) -> Dict[str, np.ndarray]:
    if not EPROP_DIR.exists():
        raise FileNotFoundError(f"Could not find E-prop directory at {EPROP_DIR}")

    np.random.seed(seed)
    sys.path.insert(0, str(EPROP_DIR))
    try:
        from input_handler import call_inputs
        from strf_handler import call_strfs
    finally:
        try:
            sys.path.remove(str(EPROP_DIR))
        except ValueError:
            pass

    p = params.as_eprop_array(input_batch_size)
    with clean_argv():
        target_dict = call_strfs(p, input_batch_size, 1)
        spks = call_inputs(1, spontaneous_fr, input_batch_size, target_dict)

    on_key = "locs_masker_None_target_0_on"
    off_key = "locs_masker_None_target_0_off"

    on_raw = spks[on_key]["stimulus_0_poisson_spks"]
    off_raw = spks[off_key]["stimulus_0_poisson_spks"]

    # Current E-prop single-cell shape is time x channel x trial x cell x batch.
    on = np.reshape(on_raw, [on_raw.shape[0], on_raw.shape[2], on_raw.shape[3], on_raw.shape[4]])
    off = np.reshape(off_raw, [off_raw.shape[0], off_raw.shape[2], off_raw.shape[3], off_raw.shape[4]])
    on = np.transpose(on, (3, 1, 2, 0))[:, :, 0, :]
    off = np.transpose(off, (3, 1, 2, 0))[:, :, 0, :]

    rate_on = spks[on_key]["stimulus_0_rate"]
    rate_off = spks[off_key]["stimulus_0_rate"]
    rate_on = np.reshape(rate_on, [rate_on.shape[0], rate_on.shape[2]])
    rate_off = np.reshape(rate_off, [rate_off.shape[0], rate_off.shape[2]])
    rate_on = np.repeat(rate_on.T[:, None, :], 10, axis=1)
    rate_off = np.repeat(rate_off.T[:, None, :], 10, axis=1)

    noise_raw = spks["noise_masker_None_target_0"]
    noise = np.transpose(noise_raw, (0, 3, 1, 2))[:, :, :, 0]

    return {
        "on_spikes": trim_or_pad(on.astype(np.float32), sim_len, axis=-1),
        "off_spikes": trim_or_pad(off.astype(np.float32), sim_len, axis=-1),
        "rate_on_hz": trim_or_pad(rate_on.astype(np.float32), sim_len, axis=-1),
        "rate_off_hz": trim_or_pad(rate_off.astype(np.float32), sim_len, axis=-1),
        "noise": trim_or_pad(noise.astype(np.float32), sim_len, axis=-1),
    }


def target_proxy_inputs(raster: np.ndarray, spontaneous_fr: float) -> Dict[str, np.ndarray]:
    trials, sim_len = raster.shape
    psth = raster.mean(axis=0, keepdims=True)
    kernel = np.ones(25, dtype=np.float32) / 25.0
    smooth = np.convolve(psth.squeeze(), kernel, mode="same")[None, :]
    on = np.repeat(smooth, trials, axis=0)[None, :, :]
    off = np.zeros_like(on)
    noise_rate = max(spontaneous_fr, 0.1) * 0.0001
    noise = np.random.binomial(1, noise_rate, size=on.shape).astype(np.float32)
    return {
        "on_spikes": on.astype(np.float32),
        "off_spikes": off.astype(np.float32),
        "rate_on_hz": on.astype(np.float32),
        "rate_off_hz": off.astype(np.float32),
        "noise": noise,
    }


def bin_raster(raster: torch.Tensor, bin_steps: int) -> torch.Tensor:
    trials, time_steps = raster.shape
    n_bins = time_steps // bin_steps
    raster = raster[:, : n_bins * bin_steps]
    return raster.reshape(trials, n_bins, bin_steps).sum(dim=(0, 2))


def sync_if_cuda(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.synchronize(device)


def reset_peak_memory_if_cuda(device: torch.device) -> None:
    if device.type == "cuda":
        torch.cuda.reset_peak_memory_stats(device)


def cuda_memory_summary(device: torch.device) -> str:
    if device.type != "cuda":
        return "CUDA memory: n/a"

    allocated = torch.cuda.memory_allocated(device) / 1024**2
    reserved = torch.cuda.memory_reserved(device) / 1024**2
    peak_allocated = torch.cuda.max_memory_allocated(device) / 1024**2
    peak_reserved = torch.cuda.max_memory_reserved(device) / 1024**2
    return (
        "CUDA memory MB "
        f"allocated={allocated:.1f}, reserved={reserved:.1f}, "
        f"peak_allocated={peak_allocated:.1f}, peak_reserved={peak_reserved:.1f}"
    )


class SingleCellLIFBackprop(nn.Module):
    dt_ms = 0.1
    e_exc = 0.0
    e_inh = -80.0

    def __init__(
        self,
        init: InitialParams,
        learn_strf_gain: bool,
        surrogate_width: float = 0.5,
    ) -> None:
        super().__init__()
        self.surrogate_width = surrogate_width

        self.raw_strf_gain = nn.Parameter(
            torch.tensor(inv_softplus(init.strf_gain), dtype=torch.float32),
            requires_grad=learn_strf_gain,
        )
        self.raw_ron_g_inc = nn.Parameter(torch.tensor(inv_softplus(init.ron_g_inc)))
        self.raw_gsyn = nn.ParameterDict(
            {
                "on_ron": nn.Parameter(torch.tensor(inv_softplus(init.on_ron_gsyn))),
                "off_ron": nn.Parameter(torch.tensor(inv_softplus(init.off_ron_gsyn))),
                "on_sonoff": nn.Parameter(torch.tensor(inv_softplus(init.on_sonoff_gsyn))),
                "off_sonoff": nn.Parameter(torch.tensor(inv_softplus(init.off_sonoff_gsyn))),
                "sonoff_ron": nn.Parameter(torch.tensor(inv_softplus(init.sonoff_ron_gsyn))),
            }
        )

    def learned_values(self) -> Dict[str, float]:
        with torch.no_grad():
            values = {
                "strf_gain": float(positive(self.raw_strf_gain).cpu()),
                "ron_g_inc": float(positive(self.raw_ron_g_inc).cpu()),
            }
            for name, raw in self.raw_gsyn.items():
                values[f"{name}_gsyn"] = float(positive(raw).cpu())
        return values

    @staticmethod
    def _filled(value: float, trials: int, device: torch.device, dtype: torch.dtype) -> torch.Tensor:
        return torch.full((trials,), value, device=device, dtype=dtype)

    def initial_state(self, trials: int, device: torch.device, dtype: torch.dtype) -> Dict[str, object]:
        z = torch.zeros(trials, device=device, dtype=dtype)
        o = torch.ones(trials, device=device, dtype=dtype)
        state: Dict[str, object] = {
            "on_v": self._filled(-65.0, trials, device, dtype),
            "off_v": self._filled(-65.0, trials, device, dtype),
            "sonoff_v": self._filled(-57.0, trials, device, dtype),
            "ron_v": self._filled(-65.0, trials, device, dtype),
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

        delay_steps = {
            "on_ron": 10,
            "off_ron": 10,
            "on_sonoff": 30,
            "off_sonoff": 30,
            "sonoff_ron": 5,
        }
        for name, steps in delay_steps.items():
            state[f"{name}_delay"] = [z.clone() for _ in range(steps)]
        return state

    @staticmethod
    def detach_state(state: Dict[str, object]) -> Dict[str, object]:
        detached: Dict[str, object] = {}
        for key, value in state.items():
            if isinstance(value, torch.Tensor):
                detached[key] = value.detach()
            elif isinstance(value, list):
                detached[key] = [item.detach() for item in value]
            else:
                detached[key] = value
        return detached

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
        if not isinstance(buffer, list) or len(buffer) == 0:
            return spike
        event = buffer.pop(0)
        buffer.append(spike)
        return event

    def _apply_synapse_event(
        self,
        state: Dict[str, object],
        name: str,
        event: torch.Tensor,
        f_p: float,
        f_f: float = 0.0,
        max_f: float = 4.0,
    ) -> None:
        x = state[f"{name}_x"]
        f = state[f"{name}_f"]
        p = state[f"{name}_p"]
        q = state[f"{name}_q"]

        state[f"{name}_x"] = x + q * event
        state[f"{name}_q"] = torch.where(event.detach() > 0, f * p, q)
        state[f"{name}_f"] = f + event * f_f * (max_f - f)
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
        trials, time_steps = on_input.shape
        if state is None:
            state = self.initial_state(trials, on_input.device, on_input.dtype)

        g = {name: positive(raw) for name, raw in self.raw_gsyn.items()}
        strf_gain = positive(self.raw_strf_gain)
        ron_g_inc = positive(self.raw_ron_g_inc)

        psth_bins = []
        trace = [] if return_trace else None
        bin_accum = torch.zeros(trials, device=on_input.device, dtype=on_input.dtype)

        for step in range(time_steps):
            on_t = on_input[:, step] * strf_gain
            off_t = off_input[:, step] * strf_gain
            noise_t = noise_input[:, step]

            self._decay_synapse(state, "on_ron", 0.7, 1.5, 30.0, 1.9481350796278847)
            self._decay_synapse(state, "off_ron", 0.7, 1.5, 30.0, 1.9481350796278847)
            self._decay_synapse(state, "on_sonoff", 0.1, 1.0, 80.0, 1.2915496650148839)
            self._decay_synapse(state, "off_sonoff", 0.1, 1.0, 80.0, 1.2915496650148839)
            self._decay_synapse(state, "sonoff_ron", 1.0, 4.5, 120.0, 1.5368523544529802)

            ron_noise_s = state["ron_noise_s"]
            ron_noise_x = state["ron_noise_x"]
            state["ron_noise_s"] = ron_noise_s + self.dt_ms * ((1.9481350796278847 * ron_noise_x - ron_noise_s) / 0.7)
            state["ron_noise_x"] = ron_noise_x + self.dt_ms * (-(ron_noise_x / 1.5) + noise_t / self.dt_ms)

            on_dv = ((-65.0 - state["on_v"]) - 200.0 * state["on_gad"] * (state["on_v"] + 80.0)
                     - 200.0 * 0.17 * on_t * (state["on_v"] - self.e_exc)) / 20.0
            off_dv = ((-65.0 - state["off_v"]) - 200.0 * state["off_gad"] * (state["off_v"] + 80.0)
                      - 200.0 * 0.17 * off_t * (state["off_v"] - self.e_exc)) / 20.0

            state["on_v"], state["on_gad"], state["on_ref"], on_spike = self._spike_neuron(
                state["on_v"], state["on_gad"], state["on_ref"], on_dv, 10.0,
                torch.as_tensor(0.0003, device=on_input.device, dtype=on_input.dtype), -47.0, -54.0, 1.0
            )
            state["off_v"], state["off_gad"], state["off_ref"], off_spike = self._spike_neuron(
                state["off_v"], state["off_gad"], state["off_ref"], off_dv, 10.0,
                torch.as_tensor(0.0003, device=on_input.device, dtype=on_input.dtype), -47.0, -54.0, 1.0
            )

            sonoff_syn = (
                g["on_sonoff"] * state["on_sonoff_s"] * (state["sonoff_v"] - self.e_exc)
                + g["off_sonoff"] * state["off_sonoff_s"] * (state["sonoff_v"] - self.e_exc)
            )
            sonoff_dv = ((-57.0 - state["sonoff_v"]) - 100.0 * state["sonoff_gad"] * (state["sonoff_v"] + 80.0)
                         - 100.0 * sonoff_syn) / 10.0
            state["sonoff_v"], state["sonoff_gad"], state["sonoff_ref"], sonoff_spike = self._spike_neuron(
                state["sonoff_v"], state["sonoff_gad"], state["sonoff_ref"], sonoff_dv, 100.0,
                torch.as_tensor(0.0, device=on_input.device, dtype=on_input.dtype), -47.0, -52.0, 0.5
            )

            ron_syn = (
                g["on_ron"] * state["on_ron_s"] * (state["ron_v"] - self.e_exc)
                + g["off_ron"] * state["off_ron_s"] * (state["ron_v"] - self.e_exc)
                + g["sonoff_ron"] * state["sonoff_ron_s"] * (state["ron_v"] - self.e_inh)
            )
            noise_current = -200.0 * 0.015 * ron_noise_s * (state["ron_v"] - self.e_exc)
            ron_dv = ((-65.0 - state["ron_v"]) - 200.0 * state["ron_gad"] * (state["ron_v"] + 80.0)
                      - 200.0 * ron_syn) / 20.0 + noise_current / 20.0
            state["ron_v"], state["ron_gad"], state["ron_ref"], ron_spike = self._spike_neuron(
                state["ron_v"], state["ron_gad"], state["ron_ref"], ron_dv, 100.0, ron_g_inc, -47.0, -54.0, 4.0
            )

            self._apply_synapse_event(state, "on_ron", self._delayed_event(state, "on_ron", on_spike), f_p=0.1)
            self._apply_synapse_event(state, "off_ron", self._delayed_event(state, "off_ron", off_spike), f_p=0.1)
            self._apply_synapse_event(state, "on_sonoff", self._delayed_event(state, "on_sonoff", on_spike), f_p=0.2)
            self._apply_synapse_event(state, "off_sonoff", self._delayed_event(state, "off_sonoff", off_spike), f_p=0.0)
            self._apply_synapse_event(state, "sonoff_ron", self._delayed_event(state, "sonoff_ron", sonoff_spike), f_p=0.5)

            bin_accum = bin_accum + ron_spike
            if return_trace:
                trace.append(ron_spike)
            if (step + 1) % bin_steps == 0:
                psth_bins.append(bin_accum.sum())
                bin_accum = torch.zeros_like(bin_accum)

        if not psth_bins:
            psth = torch.empty(0, device=on_input.device, dtype=on_input.dtype)
        else:
            psth = torch.stack(psth_bins)

        trace_tensor = torch.stack(trace, dim=1) if return_trace and trace else None
        return psth, state, trace_tensor


def select_input_arrays(
    inputs: Dict[str, np.ndarray],
    input_mode: str,
    batch_index: int,
    dt_ms: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if input_mode == "eprop-spikes" or input_mode == "target-proxy":
        on = inputs["on_spikes"][batch_index]
        off = inputs["off_spikes"][batch_index]
    elif input_mode == "eprop-rates":
        dt_s = dt_ms / 1000.0
        on = np.clip(inputs["rate_on_hz"][batch_index] * dt_s, 0.0, 1.0)
        off = np.clip(inputs["rate_off_hz"][batch_index] * dt_s, 0.0, 1.0)
    else:
        raise ValueError(f"Unknown input mode: {input_mode}")
    noise = inputs["noise"][batch_index]
    return on.astype(np.float32), off.astype(np.float32), noise.astype(np.float32)


def prepare_training_tensors(
    on: np.ndarray,
    off: np.ndarray,
    noise: np.ndarray,
    target: np.ndarray,
    max_steps: Optional[int],
    bin_steps: int,
    device: torch.device,
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    if max_steps is not None:
        on = on[:, :max_steps]
        off = off[:, :max_steps]
        noise = noise[:, :max_steps]
        target = target[:, :max_steps]

    remainder = on.shape[1] % bin_steps
    if remainder:
        on = on[:, remainder:]
        off = off[:, remainder:]
        noise = noise[:, remainder:]
        target = target[:, remainder:]

    dtype = torch.float32
    return (
        torch.as_tensor(on, dtype=dtype, device=device),
        torch.as_tensor(off, dtype=dtype, device=device),
        torch.as_tensor(noise, dtype=dtype, device=device),
        torch.as_tensor(target, dtype=dtype, device=device),
    )


def capture_output_trace(
    model: SingleCellLIFBackprop,
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
    if trace is None:
        return np.zeros((on_t.shape[0], on_t.shape[1]), dtype=np.float32)
    return trace.detach().cpu().numpy().astype(np.float32)


def train_full_bptt(
    model: SingleCellLIFBackprop,
    optimizer: torch.optim.Optimizer,
    on_t: torch.Tensor,
    off_t: torch.Tensor,
    noise_t: torch.Tensor,
    target_t: torch.Tensor,
    epochs: int,
    bin_steps: int,
    profile: bool,
) -> Tuple[list, Dict[str, object]]:
    target_bins = bin_raster(target_t, bin_steps)
    history = []
    best_record: Dict[str, object] = {
        "loss": np.inf,
        "epoch": 0,
        "output": np.zeros((target_t.shape[0], target_t.shape[1]), dtype=np.float32),
        "learned_values": model.learned_values(),
    }

    for epoch in range(epochs):
        reset_peak_memory_if_cuda(on_t.device)
        sync_if_cuda(on_t.device)
        start = time.perf_counter()

        zero_start = time.perf_counter()
        optimizer.zero_grad(set_to_none=True)
        sync_if_cuda(on_t.device)
        zero_elapsed = time.perf_counter() - zero_start

        forward_start = time.perf_counter()
        pred_bins, _, _ = model(on_t, off_t, noise_t, bin_steps)
        loss = F.mse_loss(pred_bins, target_bins)
        sync_if_cuda(on_t.device)
        forward_elapsed = time.perf_counter() - forward_start

        loss_value = float(loss.detach().cpu())
        backward_start = time.perf_counter()
        loss.backward()
        sync_if_cuda(on_t.device)
        backward_elapsed = time.perf_counter() - backward_start

        best_trace_elapsed = 0.0
        if loss_value < float(best_record["loss"]):
            best_trace_start = time.perf_counter()
            best_record = {
                "loss": loss_value,
                "epoch": epoch + 1,
                "output": capture_output_trace(model, on_t, off_t, noise_t, bin_steps),
                "learned_values": model.learned_values(),
            }
            sync_if_cuda(on_t.device)
            best_trace_elapsed = time.perf_counter() - best_trace_start

        step_start = time.perf_counter()
        optimizer.step()
        sync_if_cuda(on_t.device)
        step_elapsed = time.perf_counter() - step_start

        elapsed = time.perf_counter() - start
        history.append(loss_value)
        print(f"Epoch {epoch + 1:04d} | PSTH MSE {history[-1]:.4f} | {elapsed:.2f} s")
        if profile:
            print(
                "  profile | "
                f"zero_grad={zero_elapsed:.2f}s, forward={forward_elapsed:.2f}s, "
                f"backward={backward_elapsed:.2f}s, best_trace={best_trace_elapsed:.2f}s, "
                f"optimizer={step_elapsed:.2f}s | {cuda_memory_summary(on_t.device)}"
            )
    return history, best_record


def train_tbptt(
    model: SingleCellLIFBackprop,
    optimizer: torch.optim.Optimizer,
    on_t: torch.Tensor,
    off_t: torch.Tensor,
    noise_t: torch.Tensor,
    target_t: torch.Tensor,
    epochs: int,
    bin_steps: int,
    tbptt_steps: int,
    profile: bool,
) -> Tuple[list, Dict[str, object]]:
    if tbptt_steps % bin_steps != 0:
        raise ValueError("--tbptt-steps must be a multiple of --bin-steps")

    history = []
    usable = (on_t.shape[1] // tbptt_steps) * tbptt_steps
    on_t = on_t[:, :usable]
    off_t = off_t[:, :usable]
    noise_t = noise_t[:, :usable]
    target_t = target_t[:, :usable]
    best_record: Dict[str, object] = {
        "loss": np.inf,
        "epoch": 0,
        "output": np.zeros((target_t.shape[0], target_t.shape[1]), dtype=np.float32),
        "learned_values": model.learned_values(),
    }

    for epoch in range(epochs):
        reset_peak_memory_if_cuda(on_t.device)
        sync_if_cuda(on_t.device)
        start = time.perf_counter()
        state = None
        epoch_losses = []
        zero_elapsed = 0.0
        forward_elapsed = 0.0
        backward_elapsed = 0.0
        step_elapsed = 0.0
        for start_idx in range(0, usable, tbptt_steps):
            stop_idx = start_idx + tbptt_steps
            zero_start = time.perf_counter()
            optimizer.zero_grad(set_to_none=True)
            sync_if_cuda(on_t.device)
            zero_elapsed += time.perf_counter() - zero_start

            forward_start = time.perf_counter()
            pred_bins, state, _ = model(
                on_t[:, start_idx:stop_idx],
                off_t[:, start_idx:stop_idx],
                noise_t[:, start_idx:stop_idx],
                bin_steps,
                state=state,
            )
            target_bins = bin_raster(target_t[:, start_idx:stop_idx], bin_steps)
            loss = F.mse_loss(pred_bins, target_bins)
            sync_if_cuda(on_t.device)
            forward_elapsed += time.perf_counter() - forward_start

            backward_start = time.perf_counter()
            loss.backward()
            sync_if_cuda(on_t.device)
            backward_elapsed += time.perf_counter() - backward_start

            step_start = time.perf_counter()
            optimizer.step()
            sync_if_cuda(on_t.device)
            step_elapsed += time.perf_counter() - step_start

            epoch_losses.append(float(loss.detach().cpu()))
            state = model.detach_state(state)

        elapsed = time.perf_counter() - start
        history.append(float(np.mean(epoch_losses)))
        best_trace_elapsed = 0.0
        if history[-1] < float(best_record["loss"]):
            best_trace_start = time.perf_counter()
            best_record = {
                "loss": history[-1],
                "epoch": epoch + 1,
                "output": capture_output_trace(model, on_t, off_t, noise_t, bin_steps),
                "learned_values": model.learned_values(),
            }
            sync_if_cuda(on_t.device)
            best_trace_elapsed = time.perf_counter() - best_trace_start
        print(f"Epoch {epoch + 1:04d} | TBPTT chunk PSTH MSE {history[-1]:.4f} | {elapsed:.2f} s")
        if profile:
            print(
                "  profile | "
                f"zero_grad={zero_elapsed:.2f}s, forward={forward_elapsed:.2f}s, "
                f"backward={backward_elapsed:.2f}s, best_trace={best_trace_elapsed:.2f}s, "
                f"optimizer={step_elapsed:.2f}s | {cuda_memory_summary(on_t.device)}"
            )
    return history, best_record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-id", type=int, default=7, help="1-based cell number to fit.")
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    parser.add_argument("--sim-len", type=int, default=29801)
    parser.add_argument("--dt-ms", type=float, default=0.1)
    parser.add_argument("--bin-steps", type=int, default=200, help="200 samples = 20 ms at dt=0.1 ms.")
    parser.add_argument("--max-steps", type=int, default=None, help="Limit time samples for memory smoke tests.")
    parser.add_argument("--tbptt-steps", type=int, default=None, help="Chunk length for truncated BPTT.")
    parser.add_argument("--input-batch-size", type=int, default=1, help="Batch size used only for E-prop input generation.")
    parser.add_argument("--input-batch-index", type=int, default=0)
    parser.add_argument(
        "--input-mode",
        choices=["eprop-spikes", "eprop-rates", "target-proxy"],
        default="eprop-spikes",
        help="Use sampled E-prop spikes, continuous expected-rate tokens, or a target-derived debug proxy.",
    )
    parser.add_argument("--learn-strf-gain", action="store_true", help="Train a differentiable input/STRF gain scalar.")
    parser.add_argument("--learn-strf", dest="learn_strf_gain", action="store_true", help="Alias for --learn-strf-gain.")
    parser.add_argument("--focus-mode", choices=["peak", "first"], default="peak")
    parser.add_argument("--oliver-data-dir", type=Path, default=default_oliver_data_dir())
    parser.add_argument("--results-dir", type=Path, default=SCRIPT_DIR / "results")
    parser.add_argument("--save-trace", action="store_true", help="Save final output spike trace to the .mat file.")
    parser.add_argument("--profile", action="store_true", help="Print per-epoch timing and CUDA memory statistics.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)

    init = make_initial_params(args.seed)
    target_raster, spontaneous_fr, focus = load_cell_raster(
        args.cell_id, args.oliver_data_dir, args.dt_ms, args.sim_len, args.focus_mode
    )

    print(f"Cell {args.cell_id} | focus column {focus + 1} | spontaneous FR proxy {spontaneous_fr:.3f} Hz")
    print(f"Device: {device} | input mode: {args.input_mode} | learn STRF gain: {args.learn_strf_gain}")
    print("Initial params:", json.dumps(init.__dict__, indent=2))

    if args.input_mode == "target-proxy":
        inputs = target_proxy_inputs(target_raster, spontaneous_fr)
    else:
        inputs = generate_eprop_inputs(
            spontaneous_fr=spontaneous_fr,
            params=init,
            sim_len=args.sim_len,
            input_batch_size=args.input_batch_size,
            seed=args.seed,
        )

    on_np, off_np, noise_np = select_input_arrays(inputs, args.input_mode, args.input_batch_index, args.dt_ms)
    on_t, off_t, noise_t, target_t = prepare_training_tensors(
        on_np, off_np, noise_np, target_raster, args.max_steps, args.bin_steps, device
    )

    model = SingleCellLIFBackprop(init, learn_strf_gain=args.learn_strf_gain).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    print(f"Training tensor shape: trials={on_t.shape[0]}, time={on_t.shape[1]}")
    if args.tbptt_steps is None:
        history, best_record = train_full_bptt(
            model, optimizer, on_t, off_t, noise_t, target_t, args.epochs, args.bin_steps, args.profile
        )
    else:
        history, best_record = train_tbptt(
            model, optimizer, on_t, off_t, noise_t, target_t, args.epochs, args.bin_steps, args.tbptt_steps, args.profile
        )

    args.results_dir.mkdir(parents=True, exist_ok=True)
    result_path = args.results_dir / f"torch_backprop_cell_{args.cell_id:03d}.mat"

    final_trace = None
    if args.save_trace:
        with torch.no_grad():
            _, _, trace = model(on_t, off_t, noise_t, args.bin_steps, return_trace=True)
            final_trace = None if trace is None else trace.detach().cpu().numpy()

    best_output = np.asarray(best_record["output"], dtype=np.float32)
    best_loss = np.asarray([best_record["loss"]], dtype=np.float32)
    best_epoch = np.asarray([best_record["epoch"]], dtype=np.int32)
    best_batch_id = np.asarray([args.input_batch_index + 1], dtype=np.int32)

    result = {
        "cell_id": np.array([[args.cell_id]], dtype=np.int32),
        "focus_column": np.array([[focus + 1]], dtype=np.int32),
        "loss_history": np.asarray(history, dtype=np.float32),
        "target_raster": target_raster.astype(np.float32),
        "learned_values_json": json.dumps(model.learned_values()),
        "best_learned_values_json": json.dumps(best_record["learned_values"]),
        "input_mode": args.input_mode,
        "learn_strf_gain": np.array([[int(args.learn_strf_gain)]], dtype=np.int8),
        "max_steps": np.array([[args.max_steps or args.sim_len]], dtype=np.int32),
        "bin_steps": np.array([[args.bin_steps]], dtype=np.int32),
        "best_output": best_output,
        "best_loss": best_loss,
        "best_epoch": best_epoch,
        "best_batch_id": best_batch_id,
        "best_output_per_cell": best_output[None, :, None, :],
        "best_loss_per_cell": best_loss,
        "best_epoch_per_cell": best_epoch,
        "best_batch_id_per_cell": best_batch_id,
    }
    if final_trace is not None:
        result["final_output_trace"] = final_trace.astype(np.float32)
    savemat(result_path, result, do_compression=True)
    print("Learned values:", json.dumps(model.learned_values(), indent=2))
    print(
        f"Best output: epoch {int(best_epoch[0])}, "
        f"batch {int(best_batch_id[0])}, loss {float(best_loss[0]):.4f}"
    )
    print(f"Saved {result_path}")


if __name__ == "__main__":
    main()
