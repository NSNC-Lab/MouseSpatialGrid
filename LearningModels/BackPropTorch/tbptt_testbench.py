"""
Focused truncated-BPTT test bench for the single-cell PyTorch model.

This is intentionally smaller than benchmark_eprop_vs_backprop.py. It is for
debugging whether truncated BPTT can optimize the PSTH objective without
collapsing to silence.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

try:
    import torch
except ModuleNotFoundError as exc:
    raise SystemExit("This test bench needs PyTorch.") from exc

try:
    from scipy.io import savemat
except ModuleNotFoundError as exc:
    raise SystemExit("This test bench needs scipy.") from exc

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None

from benchmark_eprop_vs_backprop import (
    BatchedSingleCellLIFBackprop,
    bin_raster_batched,
    generate_inputs,
    load_cell_params,
    make_default_long_run_path,
    sync_torch,
)
from run_single_cell_backprop import default_oliver_data_dir, load_cell_raster


def resolve_device(name: str) -> torch.device:
    if name == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if name == "cuda" and not torch.cuda.is_available():
        raise SystemExit("CUDA was requested, but torch.cuda.is_available() is false.")
    return torch.device(name)


def capture_trace(
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


def psth_sse_by_batch(pred_bins: torch.Tensor, target_bins: torch.Tensor) -> torch.Tensor:
    return ((pred_bins - target_bins) ** 2).sum(dim=1)


def psth_from_raster(raster: np.ndarray, bin_steps: int) -> np.ndarray:
    if raster.ndim != 2:
        raise ValueError(f"Expected raster with shape trials x time, got {raster.shape}")
    n_bins = raster.shape[-1] // bin_steps
    raster = raster[:, : n_bins * bin_steps]
    return raster.reshape(raster.shape[0], n_bins, bin_steps).sum(axis=(0, 2))


def plot_raster(ax, raster: np.ndarray, color: str, dt_ms: float, title: str) -> None:
    n_trials = raster.shape[0]
    times_by_trial = [
        np.flatnonzero(raster[trial] > 0.5) * dt_ms / 1000.0 for trial in range(n_trials)
    ]
    ax.eventplot(
        times_by_trial,
        orientation="horizontal",
        lineoffsets=np.arange(1, n_trials + 1),
        linelengths=0.75,
        linewidths=0.7,
        colors=color,
    )
    ax.set_ylim(0.5, n_trials + 0.5)
    ax.invert_yaxis()
    ax.set_yticks([1, n_trials])
    ax.set_ylabel("trial", fontsize=8)
    ax.set_title(title, fontsize=9, loc="left")
    ax.tick_params(labelsize=8, length=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_epoch_top_outputs(
    path: Path,
    target_raster: np.ndarray,
    trace: np.ndarray,
    full_sse: np.ndarray,
    bin_steps: int,
    dt_ms: float,
    epoch: int,
    top_n: int,
) -> List[int]:
    if plt is None or top_n <= 0:
        return []

    top_indices = np.argsort(full_sse)[: min(top_n, trace.shape[0])]
    n_rows = len(top_indices)
    if n_rows == 0:
        return []

    target_for_plot = target_raster[:, : trace.shape[-1]]
    target_psth = psth_from_raster(target_for_plot, bin_steps)
    bin_times_s = np.arange(target_psth.size) * bin_steps * dt_ms / 1000.0

    fig, axes = plt.subplots(
        n_rows,
        3,
        figsize=(14.0, 2.1 * n_rows + 0.7),
        gridspec_kw={"width_ratios": [1.1, 1.1, 1.55]},
        squeeze=False,
    )
    fig.suptitle(
        f"Cell TBPTT epoch {epoch:04d}: top {n_rows} batches by full PSTH SSE",
        fontsize=13,
        fontweight="bold",
        y=0.995,
    )

    for row, batch_idx in enumerate(top_indices):
        sim_raster = trace[batch_idx]
        sim_psth = psth_from_raster(sim_raster, bin_steps)

        plot_raster(axes[row, 0], target_for_plot, "#111827", dt_ms, "data raster")
        plot_raster(
            axes[row, 1],
            sim_raster,
            "#2563eb",
            dt_ms,
            f"batch {batch_idx + 1} sim raster",
        )
        ax = axes[row, 2]
        ax.plot(bin_times_s, target_psth, color="#111827", lw=1.4, label="data")
        ax.plot(bin_times_s, sim_psth, color="#2563eb", lw=1.3, label="sim")
        ax.set_title(f"batch {batch_idx + 1} | SSE {full_sse[batch_idx]:.1f}", fontsize=9, loc="left")
        ax.set_ylabel("spikes / bin", fontsize=8)
        ax.tick_params(labelsize=8, length=2)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if row == 0:
            ax.legend(frameon=False, fontsize=8, loc="upper right")
        if row < n_rows - 1:
            for col in range(3):
                axes[row, col].set_xlabel("")
                axes[row, col].tick_params(labelbottom=False)
        else:
            axes[row, 0].set_xlabel("time (s)", fontsize=8)
            axes[row, 1].set_xlabel("time (s)", fontsize=8)
            ax.set_xlabel("time (s)", fontsize=8)

    fig.tight_layout(rect=(0, 0, 1, 0.975))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return [int(idx) + 1 for idx in top_indices]


def evaluate_full_psth_sse(
    model: BatchedSingleCellLIFBackprop,
    on_t: torch.Tensor,
    off_t: torch.Tensor,
    noise_t: torch.Tensor,
    target_bins: torch.Tensor,
    bin_steps: int,
) -> Tuple[np.ndarray, np.ndarray]:
    was_training = model.training
    model.eval()
    with torch.no_grad():
        pred_bins, _, trace = model(on_t, off_t, noise_t, bin_steps, return_trace=True)
        sse = psth_sse_by_batch(pred_bins, target_bins)
    if was_training:
        model.train()
    return sse.detach().cpu().numpy().astype(np.float32), trace.detach().cpu().numpy().astype(np.float32)


def train_one_epoch_tbptt(
    model: BatchedSingleCellLIFBackprop,
    optimizer: torch.optim.Optimizer,
    on_t: torch.Tensor,
    off_t: torch.Tensor,
    noise_t: torch.Tensor,
    target_t: torch.Tensor,
    bin_steps: int,
    tbptt_steps: int,
    step_mode: str,
    grad_clip: float,
) -> Dict[str, float]:
    usable = (on_t.shape[-1] // tbptt_steps) * tbptt_steps
    n_chunks = usable // tbptt_steps
    if usable <= 0 or n_chunks <= 0:
        raise ValueError("No usable TBPTT chunks. Increase --sim-len or reduce --tbptt-steps.")

    state = None
    chunk_losses: List[float] = []
    epoch_start = time.perf_counter()
    forward_s = 0.0
    backward_s = 0.0
    optimizer_s = 0.0

    if step_mode == "epoch":
        optimizer.zero_grad(set_to_none=True)

    for start_idx in range(0, usable, tbptt_steps):
        stop_idx = start_idx + tbptt_steps
        if step_mode == "chunk":
            optimizer.zero_grad(set_to_none=True)

        forward_start = time.perf_counter()
        pred_bins, state, _ = model(
            on_t[:, :, start_idx:stop_idx],
            off_t[:, :, start_idx:stop_idx],
            noise_t[:, :, start_idx:stop_idx],
            bin_steps,
            state=state,
        )
        target_bins = bin_raster_batched(target_t[:, :, start_idx:stop_idx], bin_steps)
        loss_per_batch = psth_sse_by_batch(pred_bins, target_bins)
        loss = loss_per_batch.mean()
        sync_torch(on_t.device)
        forward_s += time.perf_counter() - forward_start

        backward_start = time.perf_counter()
        loss.backward()
        sync_torch(on_t.device)
        backward_s += time.perf_counter() - backward_start

        if grad_clip > 0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip)

        if step_mode == "chunk":
            step_start = time.perf_counter()
            optimizer.step()
            sync_torch(on_t.device)
            optimizer_s += time.perf_counter() - step_start

        chunk_losses.append(float(loss.detach().cpu()))
        state = model.detach_state(state)

    if step_mode == "epoch":
        step_start = time.perf_counter()
        optimizer.step()
        sync_torch(on_t.device)
        optimizer_s += time.perf_counter() - step_start

    return {
        "chunk_sse_mean": float(np.mean(chunk_losses)),
        "chunk_sse_sum": float(np.sum(chunk_losses)),
        "forward_s": forward_s,
        "backward_s": backward_s,
        "optimizer_s": optimizer_s,
        "epoch_s": time.perf_counter() - epoch_start,
    }


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def random_eprop_params(batch_size: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    p = np.zeros((8, batch_size), dtype=np.float32)
    p[0, :] = rng.uniform(0.1, 0.2, size=batch_size).astype(np.float32)
    p[1, :] = rng.uniform(0.02, 0.03, size=batch_size).astype(np.float32)
    p[2, :] = rng.uniform(0.0, 0.00001, size=batch_size).astype(np.float32)
    p[3:5, :] = rng.uniform(0.02, 0.2, size=(2, batch_size)).astype(np.float32)
    p[5:8, :] = rng.uniform(0.0, 0.0001, size=(3, batch_size)).astype(np.float32)
    return p


def initialize_params(args: argparse.Namespace) -> Tuple[np.ndarray, str]:
    init_mode = "perturbed-longrun" if args.perturb_start else args.init_mode
    if init_mode == "random":
        return random_eprop_params(args.batch_size, args.seed), init_mode
    if init_mode == "longrun":
        return load_cell_params(args.long_run_mat, args.cell_id, args.batch_size, perturb=False), init_mode
    if init_mode == "perturbed-longrun":
        return load_cell_params(args.long_run_mat, args.cell_id, args.batch_size, perturb=True), init_mode
    raise ValueError(f"Unknown init mode: {init_mode}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cell-id", type=int, default=7)
    parser.add_argument("--batch-size", type=int, default=5)
    parser.add_argument("--epochs", type=int, default=25)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--device", choices=["auto", "cpu", "cuda"], default="auto")
    parser.add_argument("--tbptt-steps", type=int, default=100)
    parser.add_argument("--bin-steps", type=int, default=100)
    parser.add_argument("--sim-len", type=int, default=29801)
    parser.add_argument("--dt-ms", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--long-run-mat", type=Path, default=make_default_long_run_path())
    parser.add_argument("--oliver-data-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=Path(__file__).resolve().parent / "tbptt_testbench_results")
    parser.add_argument(
        "--init-mode",
        choices=["random", "longrun", "perturbed-longrun"],
        default="random",
        help="Parameter initialization source. random uses the old E-prop pinit ranges.",
    )
    parser.add_argument(
        "--step-mode",
        choices=["epoch", "chunk"],
        default="epoch",
        help="epoch accumulates chunk gradients then steps once; chunk steps after every chunk.",
    )
    parser.add_argument("--grad-clip", type=float, default=0.0, help="0 disables gradient clipping.")
    parser.add_argument("--perturb-start", action="store_true", help="Legacy alias for --init-mode perturbed-longrun.")
    parser.add_argument("--save-best-trace", action="store_true", help="Save the best full output trace in the .mat result.")
    parser.add_argument(
        "--plot-top-n",
        type=int,
        default=5,
        help="Save one raster/PSTH comparison PNG per plotted epoch for the top N batches. Use 0 to disable.",
    )
    parser.add_argument(
        "--plot-every",
        type=int,
        default=1,
        help="Save top-N raster/PSTH plots every N epochs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.tbptt_steps < args.bin_steps:
        raise SystemExit("--tbptt-steps must be >= --bin-steps for the current binned PSTH loss.")
    if args.tbptt_steps % args.bin_steps != 0:
        raise SystemExit("--tbptt-steps should be an integer multiple of --bin-steps.")
    if args.plot_every <= 0:
        raise SystemExit("--plot-every must be >= 1.")

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    device = resolve_device(args.device)

    oliver_dir = args.oliver_data_dir or default_oliver_data_dir()
    target_raster, spontaneous_fr, focus = load_cell_raster(args.cell_id, oliver_dir, args.dt_ms, args.sim_len, "peak")
    target_raster = target_raster[:, : args.sim_len].astype(np.float32)

    params, init_mode = initialize_params(args)
    inputs = generate_inputs(params, args.batch_size, spontaneous_fr, args.sim_len, args.seed)

    on_t = torch.as_tensor(inputs["on_torch"], dtype=torch.float32, device=device)
    off_t = torch.as_tensor(inputs["off_torch"], dtype=torch.float32, device=device)
    noise_t = torch.as_tensor(inputs["noise_torch"], dtype=torch.float32, device=device)
    usable = (on_t.shape[-1] // args.tbptt_steps) * args.tbptt_steps
    on_t = on_t[:, :, :usable]
    off_t = off_t[:, :, :usable]
    noise_t = noise_t[:, :, :usable]
    target = np.repeat(target_raster[None, :, :usable], args.batch_size, axis=0)
    target_t = torch.as_tensor(target, dtype=torch.float32, device=device)
    target_bins = bin_raster_batched(target_t, args.bin_steps)

    model = BatchedSingleCellLIFBackprop(params).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    args.results_dir.mkdir(parents=True, exist_ok=True)
    stem = f"tbptt_cell_{args.cell_id:03d}_batch_{args.batch_size}_{args.step_mode}_{init_mode}"
    plot_dir = args.results_dir / f"{stem}_epoch_plots"
    if args.plot_top_n > 0 and plt is None:
        print("matplotlib is not available; epoch raster/PSTH plots will be skipped.")

    print(
        f"TBPTT test bench | cell {args.cell_id} | focus {focus + 1} | "
        f"batch {args.batch_size} | device {device} | step_mode {args.step_mode} | init {init_mode}"
    )
    print(
        f"sim_len used {usable} | tbptt_steps {args.tbptt_steps} | "
        f"bin_steps {args.bin_steps} | chunks/epoch {usable // args.tbptt_steps}"
    )
    print(
        f"initial STRF gain mean {params[0].mean():.5f} | "
        f"initial STRF latency mean {params[1].mean():.5f} s"
    )

    rows: List[Dict[str, object]] = []
    initial_sse, initial_trace = evaluate_full_psth_sse(model, on_t, off_t, noise_t, target_bins, args.bin_steps)
    initial_best_idx = int(np.argmin(initial_sse))
    best = {
        "full_psth_sse_best": float(initial_sse[initial_best_idx]),
        "epoch": 0,
        "batch_id": initial_best_idx + 1,
        "trace": initial_trace[initial_best_idx].copy(),
        "params": params.copy(),
    }
    print(
        f"Epoch 0000 | full PSTH SSE mean {initial_sse.mean():.4f} | "
        f"best {initial_sse.min():.4f} batch {initial_best_idx + 1}"
    )

    for epoch in range(1, args.epochs + 1):
        train_stats = train_one_epoch_tbptt(
            model,
            optimizer,
            on_t,
            off_t,
            noise_t,
            target_t,
            args.bin_steps,
            args.tbptt_steps,
            args.step_mode,
            args.grad_clip,
        )
        full_sse, trace = evaluate_full_psth_sse(model, on_t, off_t, noise_t, target_bins, args.bin_steps)
        learned_params = model.learned_params()
        best_idx = int(np.argmin(full_sse))
        row = {
            "epoch": epoch,
            "full_psth_sse_mean": float(np.mean(full_sse)),
            "full_psth_sse_best": float(full_sse[best_idx]),
            "best_batch_id": best_idx + 1,
            "strf_gain_mean": float(np.mean(learned_params[0])),
            "strf_latency_mean_s": float(np.mean(learned_params[1])),
            "strf_gain_best_batch": float(learned_params[0, best_idx]),
            "strf_latency_best_batch_s": float(learned_params[1, best_idx]),
            "chunk_sse_mean": train_stats["chunk_sse_mean"],
            "chunk_sse_sum": train_stats["chunk_sse_sum"],
            "epoch_s": train_stats["epoch_s"],
            "forward_s": train_stats["forward_s"],
            "backward_s": train_stats["backward_s"],
            "optimizer_s": train_stats["optimizer_s"],
            "top_plot_path": "",
            "top_batch_ids": "",
        }
        if args.plot_top_n > 0 and epoch % args.plot_every == 0:
            plot_path = plot_dir / f"{stem}_epoch_{epoch:04d}_top{args.plot_top_n}.png"
            top_batch_ids = plot_epoch_top_outputs(
                plot_path,
                target_raster[:, :usable],
                trace,
                full_sse,
                args.bin_steps,
                args.dt_ms,
                epoch,
                args.plot_top_n,
            )
            if top_batch_ids:
                row["top_plot_path"] = str(plot_path)
                row["top_batch_ids"] = ",".join(str(batch_id) for batch_id in top_batch_ids)
        rows.append(row)
        if row["full_psth_sse_best"] < best["full_psth_sse_best"]:
            best = {
                "full_psth_sse_best": row["full_psth_sse_best"],
                "epoch": epoch,
                "batch_id": best_idx + 1,
                "trace": trace[best_idx].copy(),
                "params": learned_params.copy(),
            }
        print(
            f"Epoch {epoch:04d} | full PSTH SSE mean {row['full_psth_sse_mean']:.4f} | "
            f"best {row['full_psth_sse_best']:.4f} batch {row['best_batch_id']} | "
            f"STRF gain {row['strf_gain_mean']:.4f} lat {row['strf_latency_mean_s']:.5f}s | "
            f"chunk mean {row['chunk_sse_mean']:.4f} | {row['epoch_s']:.2f} s"
        )

    csv_path = args.results_dir / f"{stem}.csv"
    mat_path = args.results_dir / f"{stem}.mat"
    json_path = args.results_dir / f"{stem}_settings.json"

    write_csv(csv_path, rows)
    mat_payload: Dict[str, object] = {
        "cell_id": np.asarray([[args.cell_id]], dtype=np.int32),
        "focus_column": np.asarray([[focus + 1]], dtype=np.int32),
        "batch_size": np.asarray([[args.batch_size]], dtype=np.int32),
        "tbptt_steps": np.asarray([[args.tbptt_steps]], dtype=np.int32),
        "bin_steps": np.asarray([[args.bin_steps]], dtype=np.int32),
        "step_mode": args.step_mode,
        "init_mode": init_mode,
        "target_raster": target_raster[:, :usable],
        "initial_params": params,
        "initial_full_psth_sse": initial_sse,
        "history_full_psth_sse_mean": np.asarray([r["full_psth_sse_mean"] for r in rows], dtype=np.float32),
        "history_full_psth_sse_best": np.asarray([r["full_psth_sse_best"] for r in rows], dtype=np.float32),
        "history_best_batch_id": np.asarray([r["best_batch_id"] for r in rows], dtype=np.int32),
        "history_strf_gain_mean": np.asarray([r["strf_gain_mean"] for r in rows], dtype=np.float32),
        "history_strf_latency_mean_s": np.asarray([r["strf_latency_mean_s"] for r in rows], dtype=np.float32),
        "best_epoch": np.asarray([[best["epoch"]]], dtype=np.int32),
        "best_batch_id": np.asarray([[best["batch_id"]]], dtype=np.int32),
        "best_full_psth_sse": np.asarray([[best["full_psth_sse_best"]]], dtype=np.float32),
        "best_params": best["params"],
        "final_params": model.learned_params(),
    }
    if args.save_best_trace and best["trace"] is not None:
        mat_payload["best_output_trace"] = best["trace"]
    savemat(mat_path, mat_payload, do_compression=True)
    json_path.write_text(json.dumps(vars(args), indent=2, default=str), encoding="utf-8")

    print(f"Saved {csv_path}")
    print(f"Saved {mat_path}")


if __name__ == "__main__":
    main()
