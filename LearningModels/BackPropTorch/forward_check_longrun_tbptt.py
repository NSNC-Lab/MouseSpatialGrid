"""
Single forward-pass sanity check for the PyTorch/TBPTT model.

Loads the best final batch for one cell from LongRunResults.mat, runs the
BatchedSingleCellLIFBackprop forward pass once, and saves a raster/PSTH
comparison against the data and the original long-run E-prop output.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import torch
from scipy.io import loadmat, savemat

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmark_eprop_vs_backprop import (
    BatchedSingleCellLIFBackprop,
    bin_raster_batched,
    generate_inputs,
    make_default_long_run_path,
    sync_torch,
)
from run_single_cell_backprop import default_oliver_data_dir, load_cell_raster


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-id", type=int, default=7)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--sim-len", type=int, default=29801)
    parser.add_argument("--bin-steps", type=int, default=100)
    parser.add_argument("--dt-ms", type=float, default=0.1)
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    parser.add_argument("--long-run-mat", type=Path, default=make_default_long_run_path())
    parser.add_argument("--oliver-data-dir", type=Path, default=None)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "tbptt_forward_checks",
    )
    return parser.parse_args()


def resolve_device(name: str) -> torch.device:
    if name == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if name == "cuda" and not torch.cuda.is_available():
        raise SystemExit("CUDA was requested, but torch.cuda.is_available() is false.")
    return torch.device(name)


def load_best_longrun_cell(long_run_mat: Path, cell_id: int, sim_len: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    mat = loadmat(long_run_mat, squeeze_me=True, struct_as_record=False)
    params = np.asarray(mat["params"], dtype=np.float32)
    losses = np.asarray(mat["losses"], dtype=np.float32)
    output = np.asarray(mat["output"], dtype=np.float32)

    if params.ndim != 3:
        raise ValueError(f"Expected params with shape param x cell x batch, got {params.shape}")
    if losses.ndim != 3:
        raise ValueError(f"Expected losses with shape loss x cell x batch, got {losses.shape}")
    if output.ndim != 4:
        raise ValueError(f"Expected output with shape cell x batch x trial x time, got {output.shape}")

    best_batch = int(np.argmin(losses[1, cell_id - 1, :]))
    best_params = params[:, cell_id - 1, best_batch : best_batch + 1].copy()
    longrun_trace = output[cell_id - 1, best_batch, :, :sim_len].copy()
    best_losses = losses[:, cell_id - 1, best_batch].copy()
    return best_params, longrun_trace, best_losses, best_batch


def psth_from_raster(raster: np.ndarray, bin_steps: int) -> np.ndarray:
    n_bins = raster.shape[-1] // bin_steps
    trimmed = raster[:, : n_bins * bin_steps]
    return trimmed.reshape(trimmed.shape[0], n_bins, bin_steps).sum(axis=(0, 2))


def sse(pred_psth: np.ndarray, target_psth: np.ndarray) -> float:
    n = min(pred_psth.size, target_psth.size)
    return float(np.sum((pred_psth[:n] - target_psth[:n]) ** 2))


def plot_raster(ax, raster: np.ndarray, dt_ms: float, color: str, title: str) -> None:
    spike_times = [
        np.flatnonzero(raster[trial] > 0.5) * dt_ms / 1000.0 for trial in range(raster.shape[0])
    ]
    ax.eventplot(
        spike_times,
        orientation="horizontal",
        lineoffsets=np.arange(1, raster.shape[0] + 1),
        linelengths=0.75,
        linewidths=0.8,
        colors=color,
    )
    ax.set_ylim(0.5, raster.shape[0] + 0.5)
    ax.invert_yaxis()
    ax.set_yticks([1, raster.shape[0]])
    ax.set_ylabel("trial")
    ax.set_title(title, loc="left", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_comparison_figure(
    path: Path,
    cell_id: int,
    best_batch: int,
    target_raster: np.ndarray,
    longrun_trace: np.ndarray,
    torch_trace: np.ndarray,
    best_losses: np.ndarray,
    torch_sse: float,
    bin_steps: int,
    dt_ms: float,
    params: np.ndarray,
    device: torch.device,
) -> None:
    target_psth = psth_from_raster(target_raster, bin_steps)
    longrun_psth = psth_from_raster(longrun_trace, bin_steps)
    torch_psth = psth_from_raster(torch_trace, bin_steps)
    bin_time = np.arange(target_psth.size) * bin_steps * dt_ms / 1000.0

    fig = plt.figure(figsize=(13.5, 8.0))
    gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 1.4], hspace=0.42)
    ax_data = fig.add_subplot(gs[0])
    ax_long = fig.add_subplot(gs[1], sharex=ax_data)
    ax_torch = fig.add_subplot(gs[2], sharex=ax_data)
    ax_psth = fig.add_subplot(gs[3], sharex=ax_data)

    plot_raster(ax_data, target_raster, dt_ms, "#111827", "data raster")
    plot_raster(ax_long, longrun_trace, dt_ms, "#6b7280", f"original long-run output, batch {best_batch + 1}")
    plot_raster(ax_torch, torch_trace, dt_ms, "#2563eb", "PyTorch/TBPTT forward output")

    ax_psth.plot(bin_time, target_psth, color="#111827", lw=1.8, label="data")
    ax_psth.plot(bin_time, longrun_psth, color="#6b7280", lw=1.4, label="original long-run output")
    ax_psth.plot(bin_time, torch_psth, color="#2563eb", lw=1.4, label="PyTorch/TBPTT forward")
    ax_psth.set_ylabel("spikes / bin")
    ax_psth.set_xlabel("time (s)")
    ax_psth.spines["top"].set_visible(False)
    ax_psth.spines["right"].set_visible(False)
    ax_psth.legend(frameon=False, ncol=3, loc="upper right")

    for ax in (ax_data, ax_long, ax_torch):
        ax.tick_params(labelbottom=False)

    fig.suptitle(
        (
            f"Cell {cell_id} forward check | long-run best batch {best_batch + 1} | "
            f"device {device}"
        ),
        fontsize=14,
        fontweight="bold",
        y=0.985,
    )
    subtitle = (
        f"Long-run stored losses: L2 {best_losses[0]:.3g}, PSTH {best_losses[1]:.3g} | "
        f"PyTorch forward PSTH SSE: {torch_sse:.3g} | "
        f"params: gain {params[0,0]:.4g}, latency {params[1,0]:.4g}s"
    )
    fig.text(0.5, 0.945, subtitle, ha="center", va="center", fontsize=9, color="#374151")

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    device = resolve_device(args.device)

    params, longrun_trace, best_losses, best_batch = load_best_longrun_cell(
        args.long_run_mat, args.cell_id, args.sim_len
    )
    oliver_dir = args.oliver_data_dir or default_oliver_data_dir()
    target_raster, spontaneous_fr, focus = load_cell_raster(
        args.cell_id, oliver_dir, args.dt_ms, args.sim_len, "peak"
    )
    inputs = generate_inputs(params, 1, spontaneous_fr, args.sim_len, args.seed)

    on_t = torch.as_tensor(inputs["on_torch"], dtype=torch.float32, device=device)
    off_t = torch.as_tensor(inputs["off_torch"], dtype=torch.float32, device=device)
    noise_t = torch.as_tensor(inputs["noise_torch"], dtype=torch.float32, device=device)
    target_t = torch.as_tensor(target_raster[None, :, : on_t.shape[-1]], dtype=torch.float32, device=device)
    target_bins = bin_raster_batched(target_t, args.bin_steps)

    model = BatchedSingleCellLIFBackprop(params).to(device)
    model.eval()
    with torch.no_grad():
        pred_bins, _, trace = model(on_t, off_t, noise_t, args.bin_steps, return_trace=True)
        sync_torch(device)
        torch_sse = float(((pred_bins - target_bins) ** 2).sum().detach().cpu())
    torch_trace = trace[0].detach().cpu().numpy().astype(np.float32)

    stem = f"cell_{args.cell_id:03d}_longrun_batch_{best_batch + 1:02d}_tbptt_forward_check"
    png_path = args.output_dir / f"{stem}.png"
    mat_path = args.output_dir / f"{stem}.mat"
    json_path = args.output_dir / f"{stem}.json"

    save_comparison_figure(
        png_path,
        args.cell_id,
        best_batch,
        target_raster,
        longrun_trace,
        torch_trace,
        best_losses,
        torch_sse,
        args.bin_steps,
        args.dt_ms,
        params,
        device,
    )
    savemat(
        mat_path,
        {
            "cell_id": np.asarray([[args.cell_id]], dtype=np.int32),
            "best_batch": np.asarray([[best_batch + 1]], dtype=np.int32),
            "focus_column": np.asarray([[focus + 1]], dtype=np.int32),
            "params": params,
            "longrun_losses": best_losses,
            "torch_forward_psth_sse": np.asarray([[torch_sse]], dtype=np.float32),
            "target_raster": target_raster,
            "longrun_trace": longrun_trace,
            "torch_trace": torch_trace,
            "bin_steps": np.asarray([[args.bin_steps]], dtype=np.int32),
            "dt_ms": np.asarray([[args.dt_ms]], dtype=np.float32),
        },
        do_compression=True,
    )
    json_path.write_text(
        json.dumps(
            {
                "cell_id": args.cell_id,
                "best_batch": best_batch + 1,
                "focus_column": focus + 1,
                "long_run_mat": str(args.long_run_mat),
                "oliver_data_dir": str(oliver_dir),
                "device": str(device),
                "seed": args.seed,
                "bin_steps": args.bin_steps,
                "sim_len": args.sim_len,
                "longrun_losses": best_losses.tolist(),
                "torch_forward_psth_sse": torch_sse,
                "png_path": str(png_path),
                "mat_path": str(mat_path),
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    print(f"Saved {png_path}")
    print(f"Saved {mat_path}")
    print(f"PyTorch forward PSTH SSE: {torch_sse:.6g}")


if __name__ == "__main__":
    main()
