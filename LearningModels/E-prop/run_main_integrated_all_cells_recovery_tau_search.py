from __future__ import annotations

import csv
import importlib.util
import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat

import InitParamsAllCells
import Update_params_all_cells
import declarations
import set_options
from BuildFile import Forwards_Method_cupy_recovery_batch, calculate_loss_all_cells
from input_handler import call_inputs
from strf_handler import call_strfs


NUM_CELLS = 220
TOTAL_EPOCHS = 10
CHECKPOINT_EVERY_EPOCHS = 1

LAMBDA_VALUE = 5.0
MAX_RECOVERY_MS = 100.0
TAU_REF_LOG10_MIN = -1.0
TAU_REF_LOG10_MAX = 1.0
INITIAL_SEARCH_RADIUS_LOG10 = 1.0
MIN_SEARCH_RADIUS_LOG10 = 0.05
SEARCH_RADIUS_DECAY = 0.85

SMOOTH_COUNTS_WINDOW_STEPS = 10
FOCUS_MODE = "tuning_type"
SELECTED_TYPE = "poor+"

ISI_HIST_WEIGHT = 1.0
ISI_MEDIAN_WEIGHT = 0.25
ISI_CV_WEIGHT = 0.10
SPIKE_COUNT_WEIGHT = 0.25

ANALYSIS_DIR = Path(r"C:\Users\ipboy\Desktop\SingleChannelFigures\Figure1_MDS_Plots")
FEATURE_FIGURE_DIR = Path(r"C:\Users\ipboy\Desktop\SingleChannelFigures\Figure2_Feature_Distributions")
MODEL_OUTPUT_DIR = Path(r"C:\Users\ipboy\Documents\GitHub\ModelingEffort\Multi-Channel\Plotting\OliverDataPlotting")
ALL_DATA_PATH = ANALYSIS_DIR / "all_units_info_with_polished_criteria_modified_perf.mat"
GENERATED_SOLVER_PATH = Path(__file__).resolve().parent / "BuildFile" / "generated_solve_file_recovery_tau_search.py"


UNUSABLE = np.array(
    [
        2, 8, 11, 16, 17, 18, 20, 23, 24, 26, 34, 43, 46, 47, 48, 50, 54,
        57, 59, 60, 62, 63, 66, 70, 72, 77, 79, 81, 82, 98, 100, 105, 114,
        115, 117, 123, 125, 130, 134, 135, 137, 144, 147, 148, 151, 153,
        155, 156, 158, 159, 160, 164, 165, 166, 167, 168, 172, 173, 174,
        179, 180, 181, 182, 183, 184, 185, 187, 190, 198, 204, 206, 207,
        209, 212, 214,
    ],
    dtype=np.int32,
)
BORDERLINE_UNUSABLE = np.array(
    [
        5, 21, 22, 25, 29, 33, 37, 73, 75, 35, 55, 80, 86, 93, 94, 99,
        103, 104, 109, 119, 120, 131, 140, 152, 161, 162, 169, 170, 177,
        189, 195, 196, 197, 199, 201, 205, 211, 213,
    ],
    dtype=np.int32,
)
POOR = np.array(
    [
        3, 4, 9, 10, 27, 31, 32, 38, 39, 40, 42, 45, 49, 51, 58, 65, 67,
        69, 71, 78, 83, 85, 87, 88, 89, 95, 101, 107, 110, 111, 112, 113,
        116, 118, 121, 122, 136, 138, 141, 142, 145, 146, 150, 154, 171,
        175, 176, 193, 208, 216,
    ],
    dtype=np.int32,
)
MEDIUM = np.array(
    [
        6, 12, 13, 14, 28, 30, 36, 41, 44, 56, 61, 74, 76, 84, 90, 91,
        97, 106, 108, 126, 139, 143, 149, 157, 163, 178, 188, 191, 192,
        194, 202, 203, 215, 217, 218, 219, 220,
    ],
    dtype=np.int32,
)
GOOD = np.array([1, 15, 19, 52, 53, 64, 92, 96, 124, 127, 128, 129, 132, 186, 200, 210], dtype=np.int32)
GREAT = np.array([7, 68, 102, 133], dtype=np.int32)


def matlab_scalar_to_str(value) -> str:
    if isinstance(value, str):
        return value

    arr = np.asarray(value)
    if arr.size == 0:
        return ""

    if arr.dtype.kind in {"U", "S"}:
        if arr.ndim == 0:
            return str(arr.item())
        return "".join(arr.astype(str).ravel())

    if arr.dtype.kind == "O" and arr.size == 1:
        return matlab_scalar_to_str(arr.item())

    if arr.size == 1:
        return str(arr.item())

    return str(value)


def get_focus_index(tuning_value) -> int:
    tuning_text = matlab_scalar_to_str(tuning_value).lower()

    if "contra" in tuning_text:
        return 0
    if "45" in tuning_text:
        return 1
    if "center" in tuning_text:
        return 2
    if "ipsi" in tuning_text:
        return 3
    return 0


def choose_focus_indices(all_data, focus_mode: str, n_trials: int) -> np.ndarray:
    focus_indices = np.zeros(NUM_CELLS, dtype=np.int32)

    for cell_idx in range(NUM_CELLS):
        if focus_mode == "peak_count":
            spike_times = all_data[cell_idx].ctrl_tar1_timestamps
            counts = np.zeros(4, dtype=np.int32)
            for focus_idx in range(4):
                for trial_idx in range(n_trials):
                    times = np.asarray(spike_times[trial_idx, focus_idx]).squeeze()
                    if times.size:
                        counts[focus_idx] += np.asarray(times).size
            focus_indices[cell_idx] = int(np.argmax(counts))
        else:
            focus_indices[cell_idx] = get_focus_index(all_data[cell_idx].tuning_type)

    return focus_indices


def load_data_rasters_and_frs(all_data, opts, focus_mode: str = FOCUS_MODE):
    dt_ms = float(opts["dt"])
    time_len = int(opts["sim_len"])
    n_trials = int(opts["N_trials"])
    focus_indices = choose_focus_indices(all_data, focus_mode, n_trials)

    data_rasters = np.zeros((NUM_CELLS, n_trials, time_len), dtype=np.int8)
    spontaneous_rates = []

    for cell_idx in range(NUM_CELLS):
        spike_times = all_data[cell_idx].ctrl_tar1_timestamps[:, focus_indices[cell_idx]]
        pre_stim_spikes = []

        for trial_idx in range(n_trials):
            times = np.asarray(spike_times[trial_idx]).squeeze()
            times = np.asarray(times, dtype=np.float64).reshape(-1)
            if times.size == 0:
                continue

            stim_mask = (times >= 0) & (times < time_len * dt_ms / 1000.0)
            trial_indices = np.rint(times[stim_mask] * (1000.0 / dt_ms)).astype(np.int64)
            trial_indices = trial_indices[(trial_indices >= 0) & (trial_indices < time_len)]
            if trial_indices.size:
                data_rasters[cell_idx, trial_idx, np.unique(trial_indices)] = 1

            pre_stim_spikes.append(times[times < 0])

        pre_zeros = np.concatenate(pre_stim_spikes) if pre_stim_spikes else np.array([])
        spontaneous_rates.append(pre_zeros.size / n_trials)

    return data_rasters, spontaneous_rates, focus_indices


def moving_average(values: np.ndarray, window: int) -> np.ndarray:
    if window <= 1:
        return values.astype(np.float64, copy=True)

    values = values.astype(np.float64, copy=False)
    smoothed = np.zeros_like(values, dtype=np.float64)
    half_left = (window - 1) // 2
    half_right = window // 2

    for idx in range(values.size):
        lo = max(0, idx - half_left)
        hi = min(values.size, idx + half_right + 1)
        smoothed[idx] = np.mean(values[lo:hi])

    return smoothed


def collect_isis_steps_from_raster(raster: np.ndarray) -> np.ndarray:
    raster = np.asarray(raster)
    isis = []

    if raster.ndim == 1:
        spike_idx = np.flatnonzero(raster > 0)
        if spike_idx.size >= 2:
            isis.append(np.diff(spike_idx))
    elif raster.ndim == 2:
        for trial_idx in range(raster.shape[0]):
            spike_idx = np.flatnonzero(raster[trial_idx] > 0)
            if spike_idx.size >= 2:
                isis.append(np.diff(spike_idx))
    elif raster.ndim == 3:
        for trial_idx in range(raster.shape[0]):
            for channel_idx in range(raster.shape[1]):
                spike_idx = np.flatnonzero(raster[trial_idx, channel_idx] > 0)
                if spike_idx.size >= 2:
                    isis.append(np.diff(spike_idx))
    else:
        raise ValueError(f"Expected raster with 1, 2, or 3 dimensions, got {raster.shape}")

    if not isis:
        return np.empty((0,), dtype=np.float32)

    out = np.concatenate(isis).astype(np.float32, copy=False)
    return out[np.isfinite(out) & (out > 0)]


def estimate_vectorizing_parameters(data_rasters: np.ndarray, dt_ms: float):
    tau_abs_ms = np.full(NUM_CELLS, dt_ms, dtype=np.float32)
    rel_ref_ms = np.full(NUM_CELLS, dt_ms, dtype=np.float32)
    tau_ref_base_ms = np.full(NUM_CELLS, dt_ms, dtype=np.float32)

    bin_edges = np.arange(0, 101, 1)

    for cell_idx in range(NUM_CELLS):
        isis = collect_isis_steps_from_raster(data_rasters[cell_idx])
        if isis.size == 0:
            continue

        counts, _ = np.histogram(isis, bins=bin_edges)
        counts_smoothed = moving_average(counts, SMOOTH_COUNTS_WINDOW_STEPS)

        abs_ref = float(np.min(isis) * dt_ms)
        idx2_matlab = int(np.argmax(counts_smoothed) + 1)
        rel_ref = idx2_matlab * dt_ms - abs_ref

        tau_abs_ms[cell_idx] = max(abs_ref, 0.0)
        rel_ref_ms[cell_idx] = rel_ref
        tau_ref_base_ms[cell_idx] = max(rel_ref, dt_ms)

    return tau_abs_ms, rel_ref_ms, tau_ref_base_ms


def build_recovery_profiles(
    tau_abs_ms: np.ndarray,
    tau_ref_ms: np.ndarray,
    lambda_value: float,
    x_ms: np.ndarray,
) -> np.ndarray:
    tau_abs_ms = np.asarray(tau_abs_ms, dtype=np.float32)
    tau_ref_ms = np.asarray(tau_ref_ms, dtype=np.float32)
    x_ms = np.asarray(x_ms, dtype=np.float32)

    if tau_ref_ms.ndim == 1:
        tau_ref_ms = tau_ref_ms[:, None]

    delta = np.maximum(x_ms[None, None, :] - tau_abs_ms[:, None, None], 0.0)
    tau = np.maximum(tau_ref_ms[:, :, None], np.finfo(np.float32).eps)
    numerator = delta ** lambda_value
    denominator = numerator + tau ** lambda_value
    recovery = np.divide(numerator, denominator, out=np.zeros_like(denominator), where=denominator > 0)
    recovery = np.nan_to_num(recovery, nan=0.0, posinf=1.0, neginf=0.0)
    return np.clip(recovery, 0.0, 1.0).astype(np.float32, copy=False)


def propose_log10_scales(center_log10: np.ndarray, epoch: int, n_batch: int) -> np.ndarray:
    center_log10 = np.asarray(center_log10, dtype=np.float32).reshape(NUM_CELLS)

    if n_batch == 1:
        return np.clip(center_log10[:, None], TAU_REF_LOG10_MIN, TAU_REF_LOG10_MAX)

    if epoch == 0:
        candidates = np.linspace(TAU_REF_LOG10_MIN, TAU_REF_LOG10_MAX, n_batch, dtype=np.float32)
        return np.tile(candidates[None, :], (NUM_CELLS, 1))

    radius = max(MIN_SEARCH_RADIUS_LOG10, INITIAL_SEARCH_RADIUS_LOG10 * (SEARCH_RADIUS_DECAY ** epoch))
    offsets = np.linspace(-radius, radius, n_batch, dtype=np.float32)
    candidates = center_log10[:, None] + offsets[None, :]
    return np.clip(candidates, TAU_REF_LOG10_MIN, TAU_REF_LOG10_MAX)


def make_isi_score_bin_edges(time_len: int) -> np.ndarray:
    early_edges = np.arange(0, 501, 5, dtype=np.int32)
    late_edges = np.arange(550, time_len + 50, 50, dtype=np.int32)
    return np.unique(np.concatenate([early_edges, late_edges]))


def normalized_histogram(values: np.ndarray, bin_edges: np.ndarray) -> np.ndarray:
    counts, _ = np.histogram(values, bins=bin_edges)
    total = counts.sum()
    if total <= 0:
        return np.zeros(len(bin_edges) - 1, dtype=np.float32)
    return (counts / total).astype(np.float32)


def cv_or_nan(values: np.ndarray) -> float:
    if values.size < 2:
        return np.nan
    mean_val = float(np.mean(values))
    if mean_val <= 0:
        return np.nan
    return float(np.std(values) / mean_val)


def precompute_isi_reference(data_rasters: np.ndarray, bin_edges: np.ndarray):
    hist = np.zeros((NUM_CELLS, len(bin_edges) - 1), dtype=np.float32)
    median_steps = np.full(NUM_CELLS, np.nan, dtype=np.float32)
    cv = np.full(NUM_CELLS, np.nan, dtype=np.float32)
    n_isis = np.zeros(NUM_CELLS, dtype=np.int32)

    for cell_idx in range(NUM_CELLS):
        isis = collect_isis_steps_from_raster(data_rasters[cell_idx])
        hist[cell_idx] = normalized_histogram(isis, bin_edges)
        n_isis[cell_idx] = isis.size
        if isis.size:
            median_steps[cell_idx] = np.median(isis)
            cv[cell_idx] = cv_or_nan(isis)

    return {
        "hist": hist,
        "median_steps": median_steps,
        "cv": cv,
        "n_isis": n_isis,
    }


def safe_log_ratio(a: float, b: float) -> float:
    if not np.isfinite(a) or not np.isfinite(b):
        return 0.0
    return abs(np.log((a + 1.0) / (b + 1.0)))


def score_isi_candidates(output: np.ndarray, reference, bin_edges: np.ndarray):
    n_cells, n_batch = output.shape[:2]
    scores = np.full((n_cells, n_batch), np.inf, dtype=np.float32)
    hist_l1 = np.full((n_cells, n_batch), np.nan, dtype=np.float32)
    median_penalty = np.full((n_cells, n_batch), np.nan, dtype=np.float32)
    cv_penalty = np.full((n_cells, n_batch), np.nan, dtype=np.float32)
    count_penalty = np.full((n_cells, n_batch), np.nan, dtype=np.float32)

    for cell_idx in range(n_cells):
        ref_hist = reference["hist"][cell_idx]
        ref_median = float(reference["median_steps"][cell_idx])
        ref_cv = float(reference["cv"][cell_idx])
        ref_n = int(reference["n_isis"][cell_idx])

        for batch_idx in range(n_batch):
            isis = collect_isis_steps_from_raster(output[cell_idx, batch_idx])
            model_hist = normalized_histogram(isis, bin_edges)

            h_l1 = float(np.sum(np.abs(model_hist - ref_hist)))
            med_pen = safe_log_ratio(float(np.median(isis)) if isis.size else np.nan, ref_median)
            cv_pen = safe_log_ratio(cv_or_nan(isis), ref_cv)
            count_pen = safe_log_ratio(float(isis.size), float(ref_n))

            hist_l1[cell_idx, batch_idx] = h_l1
            median_penalty[cell_idx, batch_idx] = med_pen
            cv_penalty[cell_idx, batch_idx] = cv_pen
            count_penalty[cell_idx, batch_idx] = count_pen

            scores[cell_idx, batch_idx] = (
                ISI_HIST_WEIGHT * h_l1
                + ISI_MEDIAN_WEIGHT * med_pen
                + ISI_CV_WEIGHT * cv_pen
                + SPIKE_COUNT_WEIGHT * count_pen
            )

    return scores, {
        "hist_l1": hist_l1,
        "median_penalty": median_penalty,
        "cv_penalty": cv_penalty,
        "count_penalty": count_penalty,
    }


def build_recovery_solver(solve_file_body: str, out_path: Path = GENERATED_SOLVER_PATH):
    solve_source = (
        "# pythran export solve_run(float64[:,:,:], float64[:,:,:], float64[:,:,:], float64[:,:,:], float64[:,:]) -> Tuple[float64[:,:,:,:], float64[:,:,:]]\n"
        "import numpy as np\n"
        "import cupy as cp\n"
        "from BuildFile import calculate_loss_eprop\n"
        "def solve_run(on_input,off_input,noise_token,rate_on,rate_off,rate_on_deriv,rate_off_deriv,data,p,recovery_funcs):\n"
        + solve_file_body
    )
    out_path.write_text(solve_source, encoding="utf-8")
    print(f"{out_path.name} has been created at {out_path}.")

    module_name = "generated_solve_file_recovery_tau_search"
    spec = importlib.util.spec_from_file_location(module_name, out_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load generated recovery solver from {out_path}")

    solver = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = solver
    spec.loader.exec_module(solver)
    return solver


def selected_cell_groups(all_data):
    plot_cells = np.concatenate([POOR, MEDIUM, GOOD, GREAT]) - 1
    layer_text = [matlab_scalar_to_str(getattr(all_data[idx], "layer", "")).lower().strip() for idx in range(NUM_CELLS)]

    l23 = np.array([idx for idx in plot_cells if "2/3" in layer_text[idx]], dtype=np.int32)
    l4 = np.array([idx for idx in plot_cells if "4" in layer_text[idx]], dtype=np.int32)
    l56 = np.array([idx for idx in plot_cells if "5/6" in layer_text[idx]], dtype=np.int32)

    return [
        ("All selected", plot_cells),
        ("L2/3", l23),
        ("L4", l4),
        ("L5/6", l56),
    ]


def collect_group_isis(rasters: np.ndarray, cells: np.ndarray) -> np.ndarray:
    all_isis = []
    for cell_idx in cells:
        isis = collect_isis_steps_from_raster(rasters[cell_idx])
        if isis.size:
            all_isis.append(isis)

    if not all_isis:
        return np.empty((0,), dtype=np.float32)

    return np.concatenate(all_isis)


def add_median_marker(ax, values: np.ndarray, color, y_frac: float, x_limits):
    values = values[np.isfinite(values)]
    if values.size == 0:
        return

    med_val = float(np.median(values))
    if med_val < x_limits[0] or med_val > x_limits[1]:
        return

    y_min, y_max = ax.get_ylim()
    y_text = y_min + y_frac * (y_max - y_min)
    ax.axvline(med_val, color=color, linestyle="--", linewidth=1.2)
    ax.text(
        med_val,
        y_text,
        f" median {med_val:.1f}",
        color=color,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="center",
    )


def plot_isi_panel(ax, data_isis, model_isis, bin_edges, x_limits, dt_ms, title):
    data_color = (0.12, 0.12, 0.12)
    model_color = (0.00, 0.28, 0.85)

    if data_isis.size:
        data_weights = np.ones(data_isis.size, dtype=np.float32) / data_isis.size
        ax.hist(
            data_isis,
            bins=bin_edges,
            weights=data_weights,
            color=data_color,
            alpha=0.35,
            edgecolor="none",
        )

    model_counts, _ = np.histogram(model_isis, bins=bin_edges)
    if model_counts.sum() > 0:
        model_counts = model_counts / model_counts.sum()
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2.0
    ax.plot(bin_centers, model_counts, color=model_color, linewidth=2.0)

    add_median_marker(ax, data_isis, data_color, 0.86, x_limits)
    add_median_marker(ax, model_isis, model_color, 0.74, x_limits)

    ax.set_title(f"{title} | data n={data_isis.size}, model n={model_isis.size}", fontweight="normal")
    ax.set_xlabel(f"ISI (0.1 ms steps); 100 steps = {100 * dt_ms:.1f} ms")
    ax.set_ylabel("Probability")
    ax.set_xlim(x_limits)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_isi_group_figure(groups, data_rasters, model_rasters, bin_edges, x_limits, dt_ms, title, base_path):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    axes = axes.ravel()

    for ax, (group_name, cells) in zip(axes, groups):
        data_isis = collect_group_isis(data_rasters, cells)
        model_isis = collect_group_isis(model_rasters, cells)
        plot_isi_panel(ax, data_isis, model_isis, bin_edges, x_limits, dt_ms, group_name)

    fig.suptitle(title, fontweight="bold")
    fig.savefig(str(base_path) + ".png", dpi=300)
    fig.savefig(str(base_path) + ".pdf")
    plt.close(fig)


def mean_or_nan(values: np.ndarray) -> float:
    return float(np.mean(values)) if values.size else np.nan


def median_or_nan(values: np.ndarray) -> float:
    return float(np.median(values)) if values.size else np.nan


def save_isi_stats(groups, data_rasters, model_rasters, dt_ms, csv_path: Path):
    rows = []
    for group_name, cells in groups:
        data_isis = collect_group_isis(data_rasters, cells)
        model_isis = collect_group_isis(model_rasters, cells)

        median_data = median_or_nan(data_isis)
        median_model = median_or_nan(model_isis)
        mean_data = mean_or_nan(data_isis)
        mean_model = mean_or_nan(model_isis)

        rows.append(
            {
                "group": group_name,
                "n_data": int(data_isis.size),
                "n_model": int(model_isis.size),
                "median_data_steps": median_data,
                "median_model_steps": median_model,
                "median_data_ms": median_data * dt_ms,
                "median_model_ms": median_model * dt_ms,
                "mean_data_steps": mean_data,
                "mean_model_steps": mean_model,
                "mean_data_ms": mean_data * dt_ms,
                "mean_model_ms": mean_model * dt_ms,
                "cv_data": cv_or_nan(data_isis),
                "cv_model": cv_or_nan(model_isis),
            }
        )

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_python_isi_plots(all_data, data_rasters, best_output_per_cell, opts, timestamp: str):
    FEATURE_FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    dt_ms = float(opts["dt"])
    time_len = int(opts["sim_len"])
    model_rasters = np.asarray(best_output_per_cell)
    if model_rasters.ndim == 4 and model_rasters.shape[2] == 1:
        model_rasters = model_rasters[:, :, 0, :]

    groups = selected_cell_groups(all_data)

    full_bin_edges = np.arange(0, time_len + 50, 50, dtype=np.int32)
    zoom_bin_edges = np.arange(0, 501, 1, dtype=np.int32)

    full_base = FEATURE_FIGURE_DIR / f"ISI_distribution_full_range_python_tau_search_{SELECTED_TYPE}_{timestamp}"
    zoom_base = FEATURE_FIGURE_DIR / f"ISI_distribution_0_to_100_steps_python_tau_search_{SELECTED_TYPE}_{timestamp}"
    stats_path = FEATURE_FIGURE_DIR / f"ISI_distribution_stats_python_tau_search_{SELECTED_TYPE}_{timestamp}.csv"

    save_isi_group_figure(
        groups,
        data_rasters,
        model_rasters,
        full_bin_edges,
        (0, time_len),
        dt_ms,
        f"ISI Distribution: Full Range ({SELECTED_TYPE})",
        full_base,
    )
    save_isi_group_figure(
        groups,
        data_rasters,
        model_rasters,
        zoom_bin_edges,
        (0, 500),
        dt_ms,
        f"ISI Distribution: 0 to 500 Steps / 0 to 50 ms ({SELECTED_TYPE})",
        zoom_base,
    )
    save_isi_stats(groups, data_rasters, model_rasters, dt_ms, stats_path)

    return {
        "full_png": str(full_base) + ".png",
        "zoom_png": str(zoom_base) + ".png",
        "stats_csv": str(stats_path),
    }


def save_checkpoint(path: Path, payload: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    savemat(str(path), payload, do_compression=True)
    print(f"Saved checkpoint: {path}")


def main():
    opts = set_options.options()
    dt_ms = float(opts["dt"])
    time_len = int(opts["sim_len"])
    n_batch = int(opts["N_batch"])

    mat = loadmat(str(ALL_DATA_PATH), variable_names=["all_data"], squeeze_me=True, struct_as_record=False)
    all_data = mat["all_data"]

    data_rasters, spontaneous_rates, focus_indices = load_data_rasters_and_frs(all_data, opts)
    data = data_rasters.astype(np.float32)[:, :, :, None]

    tau_abs_ms, rel_ref_ms, tau_ref_base_ms = estimate_vectorizing_parameters(data_rasters, dt_ms)
    max_recovery_steps = round(MAX_RECOVERY_MS / dt_ms)
    recovery_x_ms = np.linspace(0, MAX_RECOVERY_MS, max_recovery_steps + 1, dtype=np.float32)

    score_bin_edges = make_isi_score_bin_edges(time_len)
    isi_reference = precompute_isi_reference(data_rasters, score_bin_edges)

    arch = declarations.Declare_Architecture(opts)
    file_body_forwards = Forwards_Method_cupy_recovery_batch.Euler_Compiler(
        arch[0],
        arch[1],
        arch[2],
        opts,
        NUM_CELLS,
    )
    generated_solve_file = build_recovery_solver(file_body_forwards)

    num_params = 8
    p, lr = InitParamsAllCells.pinit(n_batch, num_params, NUM_CELLS)
    beta1, beta2 = 0.5, 0.9995
    eps = 1e-6
    m = np.zeros((num_params, NUM_CELLS, n_batch))
    v = np.zeros((num_params, NUM_CELLS, n_batch))
    t = 0

    current_center_log10 = np.zeros(NUM_CELLS, dtype=np.float32)

    best_output_per_cell = np.zeros((NUM_CELLS, opts["N_trials"], opts["N_channels"], time_len), dtype=np.int8)
    best_isi_loss_per_cell = np.full(NUM_CELLS, np.inf, dtype=np.float32)
    best_batch_id_per_cell = np.zeros(NUM_CELLS, dtype=np.int32)
    best_epoch_per_cell = np.zeros(NUM_CELLS, dtype=np.int32)
    best_params_per_cell = np.zeros((num_params, NUM_CELLS), dtype=np.float32)
    best_tau_ref_ms_per_cell = np.full(NUM_CELLS, np.nan, dtype=np.float32)
    best_tau_ref_log10_scale_per_cell = np.full(NUM_CELLS, np.nan, dtype=np.float32)
    best_recovery_funcs_per_cell = np.zeros((NUM_CELLS, max_recovery_steps + 1), dtype=np.float32)
    best_psth_loss_at_best_isi = np.full(NUM_CELLS, np.nan, dtype=np.float32)

    tau_ref_candidates_history = np.zeros((TOTAL_EPOCHS, NUM_CELLS, n_batch), dtype=np.float32)
    tau_ref_log10_history = np.zeros((TOTAL_EPOCHS, NUM_CELLS, n_batch), dtype=np.float32)
    isi_score_history = np.zeros((TOTAL_EPOCHS, NUM_CELLS, n_batch), dtype=np.float32)
    epoch_best_batch_history = np.zeros((TOTAL_EPOCHS, NUM_CELLS), dtype=np.int32)
    epoch_best_tau_ref_history = np.zeros((TOTAL_EPOCHS, NUM_CELLS), dtype=np.float32)
    epoch_best_isi_loss_history = np.zeros((TOTAL_EPOCHS, NUM_CELLS), dtype=np.float32)

    losses_l2 = []
    losses_psth = []
    param_tracker = []

    start = time.perf_counter()
    print(
        f"Starting tau_ref ISI search with {NUM_CELLS} cells, {n_batch} tau candidates per cell, "
        f"lambda={LAMBDA_VALUE}, epochs={TOTAL_EPOCHS}."
    )

    for epoch in range(TOTAL_EPOCHS):
        tau_ref_log10_scales = propose_log10_scales(current_center_log10, epoch, n_batch)
        tau_ref_candidates_ms = tau_ref_base_ms[:, None] * (10.0 ** tau_ref_log10_scales)
        recovery_funcs = build_recovery_profiles(
            tau_abs_ms=tau_abs_ms,
            tau_ref_ms=tau_ref_candidates_ms,
            lambda_value=LAMBDA_VALUE,
            x_ms=recovery_x_ms,
        )

        tau_ref_candidates_history[epoch] = tau_ref_candidates_ms
        tau_ref_log10_history[epoch] = tau_ref_log10_scales

        target_dict = call_strfs(p, n_batch, NUM_CELLS)
        spks = call_inputs(NUM_CELLS, spontaneous_rates, n_batch, target_dict)

        on_spks = np.transpose(
            spks["locs_masker_None_target_0_on"]["stimulus_0_poisson_spks"],
            (3, 4, 2, 1, 0),
        )
        off_spks = np.transpose(
            spks["locs_masker_None_target_0_off"]["stimulus_0_poisson_spks"],
            (3, 4, 2, 1, 0),
        )
        rate_on = spks["locs_masker_None_target_0_on"]["stimulus_0_rate"]
        rate_off = spks["locs_masker_None_target_0_off"]["stimulus_0_rate"]
        rate_on_deriv = spks["locs_masker_None_target_0_on"]["stimulus_0_rate_deriv"]
        rate_off_deriv = spks["locs_masker_None_target_0_off"]["stimulus_0_rate_deriv"]
        noise = np.transpose(spks["noise_masker_None_target_0"], (0, 1, 4, 3, 2))

        output, grads, on_track, off_track, loss_holder = generated_solve_file.solve_run(
            on_spks,
            off_spks,
            noise,
            rate_on,
            rate_off,
            rate_on_deriv,
            rate_off_deriv,
            data,
            p,
            recovery_funcs,
        )

        _, loss = calculate_loss_all_cells.calculate(output, grads, data)
        losses_l2.append(np.asarray(loss[0], dtype=np.float32))
        losses_psth.append(np.asarray(loss[1], dtype=np.float32))
        param_tracker.append(np.asarray(p, dtype=np.float32).copy())

        isi_scores, isi_components = score_isi_candidates(output, isi_reference, score_bin_edges)
        isi_score_history[epoch] = isi_scores

        epoch_best_batch_idx = np.argmin(isi_scores, axis=1)
        epoch_best_isi_loss = isi_scores[np.arange(NUM_CELLS), epoch_best_batch_idx]
        epoch_best_tau_ref = tau_ref_candidates_ms[np.arange(NUM_CELLS), epoch_best_batch_idx]
        epoch_best_log10 = tau_ref_log10_scales[np.arange(NUM_CELLS), epoch_best_batch_idx]

        epoch_best_batch_history[epoch] = epoch_best_batch_idx.astype(np.int32) + 1
        epoch_best_tau_ref_history[epoch] = epoch_best_tau_ref
        epoch_best_isi_loss_history[epoch] = epoch_best_isi_loss
        current_center_log10 = epoch_best_log10.astype(np.float32, copy=False)

        improved_cells = epoch_best_isi_loss < best_isi_loss_per_cell
        if np.any(improved_cells):
            improved_idx = np.where(improved_cells)[0]
            chosen_batch_idx = epoch_best_batch_idx[improved_idx]

            best_output_per_cell[improved_idx, :, :, :] = np.asarray(
                output[improved_idx, chosen_batch_idx, :, :, :],
                dtype=np.int8,
            )
            best_isi_loss_per_cell[improved_idx] = epoch_best_isi_loss[improved_idx]
            best_batch_id_per_cell[improved_idx] = chosen_batch_idx.astype(np.int32) + 1
            best_epoch_per_cell[improved_idx] = epoch + 1
            best_params_per_cell[:, improved_idx] = p[:, improved_idx, chosen_batch_idx].astype(np.float32)
            best_tau_ref_ms_per_cell[improved_idx] = epoch_best_tau_ref[improved_idx]
            best_tau_ref_log10_scale_per_cell[improved_idx] = epoch_best_log10[improved_idx]
            best_recovery_funcs_per_cell[improved_idx] = recovery_funcs[improved_idx, chosen_batch_idx, :]

            psth_loss = np.asarray(loss[1], dtype=np.float32)
            best_psth_loss_at_best_isi[improved_idx] = psth_loss[improved_idx, chosen_batch_idx]

        print(
            f"Epoch {epoch + 1:04d}/{TOTAL_EPOCHS} | "
            f"mean PSTH loss={np.mean(loss[1]):.2f} | "
            f"median epoch ISI loss={np.nanmedian(epoch_best_isi_loss):.4f} | "
            f"median best ISI loss={np.nanmedian(best_isi_loss_per_cell):.4f} | "
            f"updated cells={int(np.sum(improved_cells))}"
        )

        grads = np.squeeze(grads)
        m, v, p, t = Update_params_all_cells.adam_update(m, v, p, t, beta1, beta2, lr, eps, grads)

        if CHECKPOINT_EVERY_EPOCHS and (epoch + 1) % CHECKPOINT_EVERY_EPOCHS == 0:
            checkpoint_payload = {
                "epoch": np.asarray(epoch + 1, dtype=np.int32),
                "params": np.asarray(p, dtype=np.float32),
                "adam_m": np.asarray(m, dtype=np.float32),
                "adam_v": np.asarray(v, dtype=np.float32),
                "adam_t": np.asarray(t, dtype=np.int32),
                "best_output_per_cell": best_output_per_cell,
                "best_isi_loss_per_cell": best_isi_loss_per_cell,
                "best_batch_id_per_cell": best_batch_id_per_cell,
                "best_epoch_per_cell": best_epoch_per_cell,
                "best_params_per_cell": best_params_per_cell,
                "best_tau_ref_ms_per_cell": best_tau_ref_ms_per_cell,
                "best_tau_ref_log10_scale_per_cell": best_tau_ref_log10_scale_per_cell,
                "best_recovery_funcs_per_cell": best_recovery_funcs_per_cell,
                "tau_abs_ms_per_cell": tau_abs_ms,
                "rel_ref_ms_per_cell": rel_ref_ms,
                "tau_ref_base_ms_per_cell": tau_ref_base_ms,
                "focus_indices_one_based": focus_indices.astype(np.int32) + 1,
                "tau_ref_candidates_history": tau_ref_candidates_history[: epoch + 1],
                "isi_score_history": isi_score_history[: epoch + 1],
                "epoch_best_batch_history": epoch_best_batch_history[: epoch + 1],
                "epoch_best_tau_ref_history": epoch_best_tau_ref_history[: epoch + 1],
                "epoch_best_isi_loss_history": epoch_best_isi_loss_history[: epoch + 1],
            }
            checkpoint_path = MODEL_OUTPUT_DIR / f"checkpoint_Eprop_All_cells_recovery_tau_search_epoch_{epoch + 1:04d}.mat"
            latest_path = MODEL_OUTPUT_DIR / "checkpoint_Eprop_All_cells_recovery_tau_search_latest.mat"
            save_checkpoint(checkpoint_path, checkpoint_payload)
            save_checkpoint(latest_path, checkpoint_payload)

    elapsed = time.perf_counter() - start
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    print(f"Finished tau_ref ISI search in {elapsed:.2f} s.")

    plot_paths = make_python_isi_plots(all_data, data_rasters, best_output_per_cell, opts, timestamp)

    final_payload = {
        "best_output_per_cell": best_output_per_cell,
        "best_isi_loss_per_cell": best_isi_loss_per_cell,
        "best_batch_id_per_cell": best_batch_id_per_cell,
        "best_epoch_per_cell": best_epoch_per_cell,
        "best_params_per_cell": best_params_per_cell,
        "best_tau_ref_ms_per_cell": best_tau_ref_ms_per_cell,
        "best_tau_ref_log10_scale_per_cell": best_tau_ref_log10_scale_per_cell,
        "best_recovery_funcs_per_cell": best_recovery_funcs_per_cell,
        "best_psth_loss_at_best_isi": best_psth_loss_at_best_isi,
        "tau_abs_ms_per_cell": tau_abs_ms,
        "rel_ref_ms_per_cell": rel_ref_ms,
        "tau_ref_base_ms_per_cell": tau_ref_base_ms,
        "lambda_value": np.asarray(LAMBDA_VALUE, dtype=np.float32),
        "recovery_x_ms": recovery_x_ms,
        "focus_indices_one_based": focus_indices.astype(np.int32) + 1,
        "losses_l2": np.asarray(losses_l2, dtype=np.float32),
        "losses_psth": np.asarray(losses_psth, dtype=np.float32),
        "params": np.asarray(param_tracker, dtype=np.float32),
        "final_params": np.asarray(p, dtype=np.float32),
        "tau_ref_candidates_history": tau_ref_candidates_history,
        "tau_ref_log10_history": tau_ref_log10_history,
        "isi_score_history": isi_score_history,
        "epoch_best_batch_history": epoch_best_batch_history,
        "epoch_best_tau_ref_history": epoch_best_tau_ref_history,
        "epoch_best_isi_loss_history": epoch_best_isi_loss_history,
        "isi_score_bin_edges": score_bin_edges,
    }

    MODEL_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = MODEL_OUTPUT_DIR / f"output_compressed_Eprop_All_cells_recovery_tau_search_{timestamp}.mat"
    savemat(str(output_path), final_payload, do_compression=True)

    best_recovery_path = ANALYSIS_DIR / f"recovery_funcs_tau_search_best_{timestamp}.mat"
    savemat(str(best_recovery_path), {"ws": best_recovery_funcs_per_cell}, do_compression=True)

    print(f"Saved final MAT output: {output_path}")
    print(f"Saved best recovery functions: {best_recovery_path}")
    print(f"Saved ISI plots/statistics: {plot_paths}")


if __name__ == "__main__":
    main()
