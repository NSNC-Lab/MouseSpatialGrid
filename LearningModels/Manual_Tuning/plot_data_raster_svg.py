import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


def _extract_data_array(mat_path: Path, data_key: str) -> np.ndarray:
    mat = loadmat(str(mat_path), squeeze_me=True, struct_as_record=False)

    if data_key in mat:
        arr = np.asarray(mat[data_key])
    else:
        candidates = []
        for key, value in mat.items():
            if key.startswith("__"):
                continue
            arr_candidate = np.asarray(value)
            if arr_candidate.ndim in (2, 3) and np.issubdtype(arr_candidate.dtype, np.number):
                candidates.append((key, arr_candidate))

        if not candidates:
            raise KeyError(
                f"Could not find key '{data_key}' and no suitable numeric 2D/3D arrays were found in {mat_path}."
            )

        arr = candidates[0][1]

    arr = np.asarray(arr)
    if arr.ndim == 3 and arr.shape[-1] == 1:
        arr = arr[:, :, 0]
    elif arr.ndim == 3 and arr.shape[0] == 1:
        arr = arr[0, :, :]

    if arr.ndim != 2:
        raise ValueError(
            f"Expected 2D trial-by-time data after squeeze, got shape {arr.shape} from {mat_path}."
        )

    return arr.astype(np.float32)


def _ensure_trials_by_time(data_2d: np.ndarray) -> np.ndarray:
    trials, timepoints = data_2d.shape
    if trials > timepoints:
        return data_2d.T
    return data_2d


def make_poster_raster(
    data_trials_time: np.ndarray,
    dt_ms: float,
    threshold: float,
    fig_width_in: float,
    fig_height_in: float,
    marker_size: float,
    line_width: float,
    title: str,
    output_svg: Path,
) -> None:
    trials, timepoints = data_trials_time.shape
    total_time_ms = float(timepoints) * dt_ms

    fig, ax = plt.subplots(figsize=(fig_width_in, fig_height_in), constrained_layout=True)

    event_times = []
    event_rows = []
    for trial_idx in range(trials):
        spike_idx = np.flatnonzero(data_trials_time[trial_idx, :] > threshold)
        if spike_idx.size == 0:
            continue
        event_times.append(spike_idx.astype(np.float32) * dt_ms)
        event_rows.append(np.full(spike_idx.shape[0], trial_idx + 1, dtype=np.int32))

    if event_times:
        x = np.concatenate(event_times)
        y = np.concatenate(event_rows)
        ax.scatter(
            x,
            y,
            s=marker_size,
            marker="|",
            linewidths=line_width,
            color="black",
        )

    ax.set_title(title, fontsize=18, pad=8)
    ax.set_xlabel("Time (ms)", fontsize=16)
    ax.set_ylabel("")

    ax.set_xlim(0, total_time_ms)
    ax.set_xticks([0, total_time_ms])
    ax.set_xticklabels(["0", f"{int(round(total_time_ms))}"])

    ax.set_yticks([])
    ax.tick_params(axis="x", which="major", labelsize=13, length=9, width=1.8, pad=5)
    ax.tick_params(axis="y", which="major", length=0)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.8)

    ax.grid(False)
    ax.set_ylim(0.5, trials + 0.5)
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    output_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_svg), format="svg", bbox_inches="tight")
    plt.close(fig)


def make_poster_psth_overlay(
    t_ms: np.ndarray,
    data_psth: np.ndarray,
    sim_psth: np.ndarray | None,
    fig_width_in: float,
    fig_height_in: float,
    title: str,
    output_svg: Path,
) -> None:
    t_ms = np.asarray(t_ms).reshape(-1)
    data_psth = np.asarray(data_psth).reshape(-1)
    sim_psth_arr = None if sim_psth is None else np.asarray(sim_psth).reshape(-1)

    fig, ax = plt.subplots(figsize=(fig_width_in, fig_height_in), constrained_layout=True)

    ax.plot(t_ms, data_psth, color="black", linewidth=1.4, label="Data")
    if sim_psth_arr is not None and sim_psth_arr.shape[0] > 0:
        ax.plot(t_ms[: sim_psth_arr.shape[0]], sim_psth_arr, color="dimgray", linewidth=1.4, label="Sim")

    total_time_ms = float(t_ms[-1]) if t_ms.shape[0] > 0 else 0.0
    ax.set_title(title, fontsize=18, pad=8)
    ax.set_xlabel("Time (ms)", fontsize=16)
    ax.set_ylabel("")

    ax.set_xlim(0, total_time_ms)
    ax.set_xticks([0, total_time_ms])
    ax.set_xticklabels(["0", f"{int(round(total_time_ms))}"])
    ax.set_yticks([])

    ax.tick_params(axis="x", which="major", labelsize=13, length=9, width=1.8, pad=5)
    ax.tick_params(axis="y", which="major", length=0)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.8)

    ax.grid(False)
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    if sim_psth_arr is not None and sim_psth_arr.shape[0] > 0:
        ax.legend(frameon=False, fontsize=11, loc="upper right")

    output_svg.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(output_svg), format="svg", bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a clean, poster-ready spike raster from data and export to SVG."
    )
    parser.add_argument(
        "--mat-file",
        type=str,
        required=True,
        help="Path to the MAT file containing trial-by-time data.",
    )
    parser.add_argument(
        "--data-key",
        type=str,
        default="picture",
        help="MAT variable key to plot (default: picture).",
    )
    parser.add_argument(
        "--output-svg",
        type=str,
        default="data_raster.svg",
        help="Output SVG path (default: data_raster.svg).",
    )
    parser.add_argument(
        "--dt-ms",
        type=float,
        default=0.1,
        help="Time step in ms per sample (default: 0.1).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.0,
        help="Values above this are treated as spikes (default: 0.0).",
    )
    parser.add_argument(
        "--fig-width-in",
        type=float,
        default=3.4,
        help="Figure width in inches (default: 3.4, good for compact poster panel).",
    )
    parser.add_argument(
        "--fig-height-in",
        type=float,
        default=1.5,
        help="Figure height in inches (default: 1.5).",
    )
    parser.add_argument(
        "--marker-size",
        type=float,
        default=6.0,
        help="Raster marker size (default: 6.0).",
    )
    parser.add_argument(
        "--line-width",
        type=float,
        default=0.8,
        help="Raster marker line width (default: 0.8).",
    )
    parser.add_argument(
        "--title",
        type=str,
        default="Data Raster",
        help="Plot title.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    mat_file = Path(args.mat_file)
    output_svg = Path(args.output_svg)

    data = _extract_data_array(mat_file, args.data_key)
    data = _ensure_trials_by_time(data)

    make_poster_raster(
        data_trials_time=data,
        dt_ms=args.dt_ms,
        threshold=args.threshold,
        fig_width_in=args.fig_width_in,
        fig_height_in=args.fig_height_in,
        marker_size=args.marker_size,
        line_width=args.line_width,
        title=args.title,
        output_svg=output_svg,
    )

    print(f"Saved SVG raster: {output_svg.resolve()}")


if __name__ == "__main__":
    main()