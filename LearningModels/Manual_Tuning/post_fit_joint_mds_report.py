from __future__ import annotations

import csv
import json
from datetime import datetime
from html import escape
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.optimize import least_squares
from scipy.stats import spearmanr
import yaml


DEFAULT_JOINT_MDS_DIR = Path(r"C:\Users\ipboy\Desktop\SingleChannelFigures\JointMDS_Space")
DEFAULT_STIMULUS_PATH = Path(
    r"C:\Users\ipboy\Documents\GitHub\MouseSpatialGrid\LearningModels\resampled-stimuli\target\200k_target1.wav"
)


def run_post_fit_joint_mds_report(
    *,
    unit_index: int,
    data: np.ndarray,
    spikes: np.ndarray,
    dt_ms: float,
    output_dir: str | Path,
    joint_mds_dir: str | Path = DEFAULT_JOINT_MDS_DIR,
    stimulus_path: str | Path = DEFAULT_STIMULUS_PATH,
    condition_label: str = "contra",
    manual_param_matrix: np.ndarray | None = None,
    manual_param_names: List[str] | None = None,
    config_path: str | Path | None = None,
    extra_summary_fields: Dict[str, object] | None = None,
) -> Dict[str, object]:
    output_root = Path(output_dir)
    output_dir = output_root / f"unit_{int(unit_index):03d}"
    output_dir.mkdir(parents=True, exist_ok=True)

    joint_mds_dir = Path(joint_mds_dir)
    stimulus_path = Path(stimulus_path)

    data_trials_time = coerce_trial_time(data)
    model_batches = coerce_model_batches(spikes)
    model_batches = model_batches[:, : data_trials_time.shape[0], : data_trials_time.shape[1]]
    data_trials_time = data_trials_time[: model_batches.shape[1], : model_batches.shape[2]]

    batch_scores = []
    for batch_idx in range(model_batches.shape[0]):
        batch_metrics = compute_fit_metrics(
            model_trial_time=model_batches[batch_idx],
            data_trial_time=data_trials_time,
            dt_ms=dt_ms,
            stimulus_path=stimulus_path,
        )
        batch_scores.append(
            {
                "batch_idx": batch_idx,
                "psth_20ms_sse": batch_metrics["psth_20ms_sse"],
                "rcorr": batch_metrics["rcorr"],
                "psth_20ms_pearson_r": batch_metrics["psth_20ms_pearson_r"],
                "noise_corrected_r_10ms": batch_metrics["noise_corrected_r_10ms"],
            }
        )

    best_batch_idx = int(np.argmin([row["psth_20ms_sse"] for row in batch_scores]))
    best_model_trial_time = model_batches[best_batch_idx]
    metrics = compute_fit_metrics(
        model_trial_time=best_model_trial_time,
        data_trial_time=data_trials_time,
        dt_ms=dt_ms,
        stimulus_path=stimulus_path,
    )

    batch_table_path = output_dir / "post_fit_batch_ranking.csv"
    write_dict_rows_csv(batch_table_path, batch_scores)

    metrics_table_path = output_dir / "post_fit_metrics.csv"
    write_metric_csv(metrics_table_path, metrics)

    reference = load_reference_joint_mds(joint_mds_dir)
    unit_zero = int(unit_index) - 1
    if unit_zero < 0 or unit_zero >= int(reference["n_cells"]):
        raise IndexError(f"Unit index {unit_index} is outside the reference Joint MDS range 1..{reference['n_cells']}.")

    model_psth_20ms = summed_psth(best_model_trial_time, ms_to_samples(20.0, dt_ms))
    manual_fit_xy = project_point_to_reference_space(model_psth_20ms, reference["psths_all"], reference["y_ref"])
    reference_data_xy = np.asarray(reference["y_ref"][unit_zero], dtype=np.float64)
    reference_model_xy = np.asarray(reference["y_ref"][reference["n_cells"] + unit_zero], dtype=np.float64)
    reference_batch_idx = int(reference["best_batch_idx"][unit_zero])

    projection_payload = {
        "reference_y": reference["y_ref"],
        "n_cells": reference["n_cells"],
        "reference_data_xy": reference_data_xy,
        "reference_model_xy": reference_model_xy,
        "manual_fit_xy": manual_fit_xy,
        "unit_index": unit_index,
        "best_batch_idx": best_batch_idx,
        "reference_best_batch_idx": reference_batch_idx,
    }

    figure_path = output_dir / "post_fit_joint_mds_report.svg"
    save_report_figure(figure_path, metrics, projection_payload)

    manual_distance = float(np.linalg.norm(manual_fit_xy - reference_data_xy))
    original_model_distance = float(np.linalg.norm(reference_model_xy - reference_data_xy))
    distance_improvement = float(original_model_distance - manual_distance)

    summary_row = {
        "run_timestamp": datetime.now().isoformat(timespec="seconds"),
        "unit_index": int(unit_index),
        "unit_folder": str(output_dir),
        "condition_label": condition_label,
        "manual_best_batch_idx": int(best_batch_idx),
        "reference_best_batch_idx": int(reference_batch_idx),
        "mds_manual_fit_distance_to_data": manual_distance,
        "mds_original_model_distance_to_data": original_model_distance,
        "mds_distance_improvement": distance_improvement,
        "reference_data_mds_x": float(reference_data_xy[0]),
        "reference_data_mds_y": float(reference_data_xy[1]),
        "reference_model_mds_x": float(reference_model_xy[0]),
        "reference_model_mds_y": float(reference_model_xy[1]),
        "manual_fit_mds_x": float(manual_fit_xy[0]),
        "manual_fit_mds_y": float(manual_fit_xy[1]),
    }
    summary_row.update(extract_best_batch_params(manual_param_matrix, best_batch_idx, manual_param_names))
    summary_row.update(load_config_parameter_summary(config_path))
    if extra_summary_fields:
        summary_row.update(extra_summary_fields)
    summary_row.update(metrics)

    cumulative_csv_path = output_root / "manual_tuning_post_fit_summary.csv"
    cumulative_html_path = output_root / "manual_tuning_post_fit_summary.html"
    update_summary_outputs(cumulative_csv_path, cumulative_html_path, summary_row)

    unit_history_csv_path = output_dir / "post_fit_unit_history.csv"
    unit_history_html_path = output_dir / "post_fit_unit_history.html"
    update_summary_outputs(unit_history_csv_path, unit_history_html_path, summary_row)

    manifest_path = output_dir / "post_fit_run_manifest.json"
    write_run_manifest(
        manifest_path,
        summary_row=summary_row,
        metrics_csv=metrics_table_path,
        batch_csv=batch_table_path,
        figure_path=figure_path,
        projection_npz=output_dir / "post_fit_projection_data.npz",
    )

    np.savez_compressed(
        output_dir / "post_fit_projection_data.npz",
        reference_data_xy=reference_data_xy,
        reference_model_xy=reference_model_xy,
        manual_fit_xy=manual_fit_xy,
        best_batch_idx=best_batch_idx,
        reference_best_batch_idx=reference_batch_idx,
        reference_y=reference["y_ref"],
    )

    print(f"Post-fit report saved to: {output_dir}")
    print(f"Best batch by 20 ms PSTH SSE: {best_batch_idx}")
    print(f"Joint MDS report figure: {figure_path}")
    print(f"Updated unit history CSV: {unit_history_csv_path}")
    print(f"Updated unit history HTML: {unit_history_html_path}")
    print(f"Updated cumulative summary CSV: {cumulative_csv_path}")
    print(f"Updated cumulative summary HTML: {cumulative_html_path}")

    return {
        "output_dir": output_dir,
        "best_batch_idx": best_batch_idx,
        "metrics": metrics,
        "projection_figure": figure_path,
        "metrics_csv": metrics_table_path,
        "batch_csv": batch_table_path,
        "unit_history_csv": unit_history_csv_path,
        "unit_history_html": unit_history_html_path,
        "manifest_json": manifest_path,
        "summary_csv": cumulative_csv_path,
        "summary_html": cumulative_html_path,
    }


def coerce_trial_time(data: np.ndarray) -> np.ndarray:
    arr = np.asarray(data, dtype=np.float32).squeeze()
    if arr.ndim == 1:
        arr = arr[None, :]
    elif arr.ndim == 3:
        if arr.shape[-1] == 1:
            arr = arr[..., 0]
        else:
            arr = arr.sum(axis=1)
    if arr.shape[0] > arr.shape[1] and arr.shape[1] <= 32:
        arr = arr.T
    return np.asarray(arr, dtype=np.float32)


def coerce_model_batches(spikes: np.ndarray) -> np.ndarray:
    arr = np.asarray(spikes, dtype=np.float32)
    if arr.ndim == 4:
        return arr.sum(axis=2)
    if arr.ndim == 3:
        if arr.shape[0] <= 64 and arr.shape[1] <= 64 and arr.shape[2] > 100:
            return arr
        return arr[None, ...]
    if arr.ndim == 2:
        return arr[None, ...]
    raise ValueError(f"Unsupported spike array shape: {arr.shape}")


def ms_to_samples(ms: float, dt_ms: float) -> int:
    width = int(round(ms / dt_ms))
    return max(width, 1)


def bin_trial_time(trial_time: np.ndarray, bin_width_samples: int) -> np.ndarray:
    n_time = trial_time.shape[1]
    usable = (n_time // bin_width_samples) * bin_width_samples
    if usable <= 0:
        return np.zeros((trial_time.shape[0], 0), dtype=np.float32)
    trimmed = trial_time[:, :usable]
    return trimmed.reshape(trial_time.shape[0], -1, bin_width_samples).sum(axis=2)


def summed_psth(trial_time: np.ndarray, bin_width_samples: int) -> np.ndarray:
    return bin_trial_time(trial_time, bin_width_samples).sum(axis=0)


def mean_rate_hz(trial_time: np.ndarray, dt_ms: float) -> float:
    duration_sec = trial_time.shape[1] * dt_ms / 1000.0
    if duration_sec <= 0:
        return np.nan
    return float(trial_time.sum(axis=1).mean() / duration_sec)


def safe_pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    y = np.asarray(y, dtype=np.float64).reshape(-1)
    if x.size < 2 or y.size < 2:
        return np.nan
    if np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def safe_cosine(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64).reshape(-1)
    y = np.asarray(y, dtype=np.float64).reshape(-1)
    denom = np.linalg.norm(x) * np.linalg.norm(y)
    if denom <= 0:
        return np.nan
    return float(np.dot(x, y) / denom)


def lifetime_sparseness(response: np.ndarray) -> float:
    response = np.asarray(response, dtype=np.float64).reshape(-1)
    response = response[np.isfinite(response)]
    n_bins = response.size
    if n_bins <= 1:
        return np.nan
    if response.size == 0:
        return np.nan
    response = np.clip(response, a_min=0.0, a_max=None)
    normfact = float(np.sum(response * response))
    if normfact <= 0.0:
        return np.nan
    numerator = n_bins - (float(np.sum(response)) ** 2) / normfact
    return float(numerator / (n_bins - 1))


def extract_best_batch_params(
    manual_param_matrix: np.ndarray | None,
    best_batch_idx: int,
    manual_param_names: List[str] | None,
) -> Dict[str, object]:
    if manual_param_matrix is None:
        return {}

    params = np.asarray(manual_param_matrix, dtype=np.float64)
    if params.ndim != 2:
        return {}
    if best_batch_idx < 0 or best_batch_idx >= params.shape[1]:
        return {}

    n_params = params.shape[0]
    if manual_param_names is None or len(manual_param_names) != n_params:
        manual_param_names = [f"manual_param_{idx}" for idx in range(n_params)]

    summary: Dict[str, object] = {}
    for idx, name in enumerate(manual_param_names):
        safe_name = sanitize_summary_key(name)
        summary[f"manual_{safe_name}"] = float(params[idx, best_batch_idx])
    return summary


def load_config_parameter_summary(config_path: str | Path | None) -> Dict[str, object]:
    if config_path is None:
        return {}

    path = Path(config_path)
    if not path.exists():
        return {"config_path": str(path)}

    with path.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}

    strf_config = config.get("strf_config", {}) or {}
    param_h = strf_config.get("paramH", {}) or {}

    summary: Dict[str, object] = {"config_path": str(path)}
    if "strfGain" in strf_config:
        summary["config_strfGain"] = float(strf_config["strfGain"])
    if "alpha" in param_h:
        summary["config_paramH_alpha"] = float(param_h["alpha"])
    return summary


def sanitize_summary_key(name: str) -> str:
    sanitized = []
    for char in str(name):
        if char.isalnum():
            sanitized.append(char)
        else:
            sanitized.append("_")
    compact = "".join(sanitized).strip("_")
    while "__" in compact:
        compact = compact.replace("__", "_")
    return compact


CANONICAL_PARAMETER_KEYS = (
    "On_ROn_gSYN",
    "Off_ROn_gSYN",
    "On_SOnOff_gSYN",
    "Off_SOnOff_gSYN",
    "SOnOff_ROn_gSYN",
    "ROn_tau_ad",
    "ROn_g_inc",
)


def update_summary_outputs(csv_path: Path, html_path: Path, new_row: Dict[str, object]) -> None:
    rows = read_csv_rows(csv_path)
    rows.append(normalize_summary_row(dict(new_row)))

    headers: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in headers:
                headers.append(key)

    write_csv_rows(csv_path, headers, rows)
    write_summary_html(html_path, headers, rows)


def read_csv_rows(path: Path) -> List[Dict[str, object]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return [normalize_summary_row(dict(row)) for row in reader]


def normalize_summary_row(row: Dict[str, object]) -> Dict[str, object]:
    normalized: Dict[str, object] = {}
    drop_keys = legacy_parameter_keys()

    for key, value in row.items():
        if key in drop_keys:
            continue
        normalized[key] = value

    for base_key in CANONICAL_PARAMETER_KEYS:
        selected = select_parameter_value(row, base_key)
        if selected is not None:
            normalized[base_key] = selected

    return normalized


def legacy_parameter_keys() -> set[str]:
    prefixes = ("generated_", "declared_", "architecture_")
    keys = set()
    for base_key in CANONICAL_PARAMETER_KEYS:
        for prefix in prefixes:
            keys.add(f"{prefix}{base_key}")
    keys.update(
        {
            "manual_On_ROn_gSYN",
            "manual_Off_ROn_gSYN",
            "manual_On_SOnOff_gSYN",
            "manual_Off_SOnOff_gSYN",
            "manual_SOnOff_ROn_gSYN",
        }
    )
    return keys


def select_parameter_value(row: Dict[str, object], base_key: str) -> object | None:
    candidate_keys = (
        base_key,
        f"generated_{base_key}",
        f"declared_{base_key}",
        f"architecture_{base_key}",
        f"manual_{base_key}",
    )
    for key in candidate_keys:
        value = row.get(key, None)
        if has_summary_value(value):
            return value
    return None


def has_summary_value(value: object) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return value.strip() != ""
    return True


def write_csv_rows(path: Path, headers: List[str], rows: List[Dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            normalized = {header: stringify_csv_value(row.get(header, "")) for header in headers}
            writer.writerow(normalized)


def stringify_csv_value(value: object) -> str:
    if isinstance(value, (float, np.floating)):
        if np.isnan(value):
            return "nan"
        return f"{float(value):.12g}"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return str(value)


def write_summary_html(path: Path, headers: List[str], rows: List[Dict[str, object]]) -> None:
    display_rows = list(reversed(rows))
    html_lines = [
        "<!doctype html>",
        "<html lang='en'>",
        "<head>",
        "<meta charset='utf-8'>",
        "<title>Manual Tuning Post-Fit Summary</title>",
        "<style>",
        "body { font-family: Arial, sans-serif; margin: 20px; background: #fafafa; color: #1f1f1f; }",
        "h1 { margin-bottom: 0.25rem; }",
        "p { margin-top: 0; color: #555; }",
        ".table-wrap { overflow-x: auto; border: 1px solid #d9d9d9; background: white; }",
        "table { border-collapse: collapse; width: max-content; min-width: 100%; }",
        "thead th { position: sticky; top: 0; background: #f0f3f8; z-index: 1; }",
        "th, td { border: 1px solid #d9d9d9; padding: 6px 8px; font-size: 12px; white-space: nowrap; }",
        "tbody tr:nth-child(even) { background: #fbfbfb; }",
        "td.num { text-align: right; font-variant-numeric: tabular-nums; }",
        "</style>",
        "</head>",
        "<body>",
        "<h1>Manual Tuning Post-Fit Summary</h1>",
        f"<p>{len(display_rows)} runs logged. Newest run is shown first.</p>",
        "<div class='table-wrap'>",
        "<table>",
        "<thead><tr>",
    ]
    html_lines.extend(f"<th>{escape(header)}</th>" for header in headers)
    html_lines.extend(["</tr></thead>", "<tbody>"])

    for row in display_rows:
        html_lines.append("<tr>")
        for header in headers:
            value = stringify_csv_value(row.get(header, ""))
            css_class = "num" if looks_numeric(value) else ""
            html_lines.append(f"<td class='{css_class}'>{escape(value)}</td>")
        html_lines.append("</tr>")

    html_lines.extend([
        "</tbody>",
        "</table>",
        "</div>",
        "</body>",
        "</html>",
    ])
    path.write_text("\n".join(html_lines) + "\n", encoding="utf-8")


def looks_numeric(value: str) -> bool:
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


def write_run_manifest(
    path: Path,
    *,
    summary_row: Dict[str, object],
    metrics_csv: Path,
    batch_csv: Path,
    figure_path: Path,
    projection_npz: Path,
) -> None:
    payload = {
        "summary_row": {key: stringify_csv_value(value) for key, value in summary_row.items()},
        "paths": {
            "metrics_csv": str(metrics_csv),
            "batch_csv": str(batch_csv),
            "figure_path": str(figure_path),
            "projection_npz": str(projection_npz),
        },
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def split_half_corr(trial_bins: np.ndarray) -> float:
    if trial_bins.shape[0] < 2:
        return np.nan
    half_a = trial_bins[::2]
    half_b = trial_bins[1::2]
    if half_b.shape[0] == 0:
        midpoint = trial_bins.shape[0] // 2
        half_a = trial_bins[:midpoint]
        half_b = trial_bins[midpoint:]
    if half_a.shape[0] == 0 or half_b.shape[0] == 0:
        return np.nan
    return safe_pearson(half_a.mean(axis=0), half_b.mean(axis=0))


def noise_corrected_r(model_trial_bins: np.ndarray, data_trial_bins: np.ndarray) -> float:
    prediction = model_trial_bins.mean(axis=0)
    trial_pred_rs = [safe_pearson(prediction, trial) for trial in data_trial_bins]

    pair_rs = []
    for first_idx in range(data_trial_bins.shape[0] - 1):
        for second_idx in range(first_idx + 1, data_trial_bins.shape[0]):
            pair_rs.append(safe_pearson(data_trial_bins[first_idx], data_trial_bins[second_idx]))

    trial_pred_rs = np.asarray(trial_pred_rs, dtype=np.float64)
    pair_rs = np.asarray(pair_rs, dtype=np.float64)

    if not np.isfinite(trial_pred_rs).any() or not np.isfinite(pair_rs).any():
        return np.nan

    trial_pred_mean = float(np.nanmean(trial_pred_rs))
    pair_mean = float(np.nanmean(pair_rs))
    pair_mean = max(pair_mean, 0.01)
    return float(trial_pred_mean / np.sqrt(pair_mean))


def get_stimulus_event_info(stimulus_path: Path, response_bin_ms: float, n_bins: int) -> Dict[str, np.ndarray]:
    from scipy.io.wavfile import read as wavread

    stimulus_fs, waveform = wavread(str(stimulus_path))
    waveform = np.asarray(waveform, dtype=np.float64)
    if waveform.ndim > 1:
        waveform = waveform.mean(axis=1)

    envelope = np.abs(waveform)
    smooth_window = max(1, int(round(0.005 * stimulus_fs)))
    envelope = moving_average(envelope, smooth_window)

    samples_per_bin = max(1, int(round((response_bin_ms / 1000.0) * stimulus_fs)))
    usable = min(envelope.shape[0], samples_per_bin * n_bins)
    envelope = envelope[:usable]
    if envelope.shape[0] < samples_per_bin * n_bins:
        envelope = np.pad(envelope, (0, samples_per_bin * n_bins - envelope.shape[0]))

    envelope_bins = envelope.reshape(n_bins, samples_per_bin).mean(axis=1)
    envelope_bins = moving_average(envelope_bins, 3)
    if envelope_bins.max() > 0:
        envelope_bins = envelope_bins / envelope_bins.max()

    threshold = 0.05
    active_mask = envelope_bins > threshold
    start_idx = int(np.argmax(active_mask)) if np.any(active_mask) else 0
    end_idx = int(np.where(active_mask)[0][-1]) if np.any(active_mask) else n_bins - 1
    first_peak_idx = find_first_local_max(envelope_bins, start_idx)
    first_trough_idx = find_first_local_min(envelope_bins, first_peak_idx)
    all_onset_bins, all_offset_bins = build_rise_fall_masks(envelope_bins, start_idx, end_idx)

    return {
        "first_onset_bins": np.arange(start_idx, max(first_peak_idx, start_idx) + 1, dtype=int),
        "first_offset_bins": np.arange(first_peak_idx, max(first_trough_idx, first_peak_idx) + 1, dtype=int),
        "onset_bins": all_onset_bins,
        "offset_bins": all_offset_bins,
    }


def moving_average(x: np.ndarray, width: int) -> np.ndarray:
    if width <= 1:
        return np.asarray(x, dtype=np.float64)
    kernel = np.ones(width, dtype=np.float64) / float(width)
    return np.convolve(x, kernel, mode="same")


def find_first_local_max(values: np.ndarray, start_idx: int) -> int:
    for idx in range(start_idx + 1, max(start_idx + 1, values.size - 1)):
        if values[idx] >= values[idx - 1] and values[idx] >= values[idx + 1] and np.any(values[start_idx : idx + 1] > values[start_idx]):
            return idx
    return int(start_idx + np.argmax(values[start_idx:]))


def find_first_local_min(values: np.ndarray, start_idx: int) -> int:
    for idx in range(start_idx + 1, max(start_idx + 1, values.size - 1)):
        if values[idx] <= values[idx - 1] and values[idx] <= values[idx + 1] and np.any(values[start_idx : idx + 1] < values[start_idx]):
            return idx
    return int(start_idx + np.argmin(values[start_idx:]))


def build_rise_fall_masks(envelope_bins: np.ndarray, start_idx: int, end_idx: int) -> Tuple[np.ndarray, np.ndarray]:
    if end_idx <= start_idx:
        return np.array([start_idx], dtype=int), np.array([], dtype=int)

    delta = np.diff(envelope_bins[start_idx : end_idx + 1])
    if delta.size == 0:
        return np.array([start_idx], dtype=int), np.array([], dtype=int)

    slope_threshold = max(1e-4, 0.01 * np.max(np.abs(delta)))
    state = np.zeros_like(delta, dtype=int)
    state[delta > slope_threshold] = 1
    state[delta < -slope_threshold] = -1
    state = fill_local_state(state)

    if not np.any(state):
        return np.arange(start_idx, end_idx + 1, dtype=int), np.array([], dtype=int)

    bin_state = np.zeros(end_idx - start_idx + 1, dtype=int)
    bin_state[:-1] = state
    bin_state[-1] = state[-1]

    onset_bins = (start_idx + np.flatnonzero(bin_state > 0)).astype(int)
    offset_bins = (start_idx + np.flatnonzero(bin_state < 0)).astype(int)
    return onset_bins, offset_bins


def fill_local_state(state: np.ndarray) -> np.ndarray:
    state = state.copy()
    nonzero = np.flatnonzero(state)
    if nonzero.size == 0:
        return state
    first = nonzero[0]
    state[:first] = state[first]
    for idx in range(first + 1, state.size):
        if state[idx] == 0:
            state[idx] = state[idx - 1]
    return state


def compute_response_features(trial_bins: np.ndarray, bin_ms: float, stimulus_info: Dict[str, np.ndarray]) -> Dict[str, float]:
    bin_sec = bin_ms / 1000.0
    mean_counts = trial_bins.mean(axis=0)
    mean_rate = mean_counts / bin_sec
    centers_ms = (np.arange(mean_rate.size) + 0.5) * bin_ms

    onset_window = clip_indices(stimulus_info["first_onset_bins"], mean_rate.size)
    offset_window = clip_indices(stimulus_info["first_offset_bins"], mean_rate.size)
    peak_rate, peak_idx = first_event_peak(mean_rate, onset_window, offset_window)

    if np.isnan(peak_rate) or peak_idx is None:
        peak_rate = np.nan
        peak_latency_ms = np.nan
    else:
        peak_latency_ms = centers_ms[peak_idx]

    onset_resp = window_responsiveness_index(mean_rate, clip_indices(stimulus_info["onset_bins"], mean_rate.size))
    offset_resp = window_responsiveness_index(mean_rate, clip_indices(stimulus_info["offset_bins"], mean_rate.size))

    return {
        "peak_rate_hz": float(peak_rate),
        "peak_latency_ms": float(peak_latency_ms),
        "onset_responsiveness": float(onset_resp),
        "offset_responsiveness": float(offset_resp),
    }


def clip_indices(indices: np.ndarray, max_len: int) -> np.ndarray:
    indices = np.asarray(indices, dtype=int)
    indices = indices[(indices >= 0) & (indices < max_len)]
    return np.unique(indices)


def first_event_peak(mean_rate: np.ndarray, onset_window: np.ndarray, offset_window: np.ndarray) -> Tuple[float, int | None]:
    onset_peak = np.nan if onset_window.size == 0 else float(np.max(mean_rate[onset_window]))
    offset_peak = np.nan if offset_window.size == 0 else float(np.max(mean_rate[offset_window]))

    if np.isfinite(onset_peak) and (not np.isfinite(offset_peak) or onset_peak >= offset_peak):
        return onset_peak, int(onset_window[np.argmax(mean_rate[onset_window])])
    if np.isfinite(offset_peak):
        return offset_peak, int(offset_window[np.argmax(mean_rate[offset_window])])
    return np.nan, None


def window_responsiveness_index(mean_rate: np.ndarray, focus_window: np.ndarray) -> float:
    if focus_window.size == 0:
        return np.nan
    focus_mask = np.zeros(mean_rate.size, dtype=bool)
    focus_mask[focus_window] = True
    background_mask = ~focus_mask
    if not np.any(background_mask):
        return np.nan
    focus_rate = float(np.mean(mean_rate[focus_mask]))
    background_rate = float(np.mean(mean_rate[background_mask]))
    return float((focus_rate - background_rate) / (focus_rate + background_rate + np.finfo(float).eps))


def compute_fit_metrics(
    *,
    model_trial_time: np.ndarray,
    data_trial_time: np.ndarray,
    dt_ms: float,
    stimulus_path: Path,
) -> Dict[str, float]:
    metrics: Dict[str, float] = {}

    model_trial_time = np.asarray(model_trial_time, dtype=np.float32)
    data_trial_time = np.asarray(data_trial_time, dtype=np.float32)
    diff_raw = model_trial_time - data_trial_time
    metrics["raw_raster_mse"] = float(np.mean(diff_raw * diff_raw))
    metrics["raw_raster_rmse"] = float(np.sqrt(metrics["raw_raster_mse"]))
    metrics["raw_raster_mae"] = float(np.mean(np.abs(diff_raw)))

    binned_cache: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
    for bin_ms in (5.0, 10.0, 20.0, 50.0):
        width = ms_to_samples(bin_ms, dt_ms)
        model_bins = bin_trial_time(model_trial_time, width)
        data_bins = bin_trial_time(data_trial_time, width)
        binned_cache[int(bin_ms)] = (model_bins, data_bins)
        model_psth = model_bins.sum(axis=0)
        data_psth = data_bins.sum(axis=0)
        diff = model_psth - data_psth
        metrics[f"psth_{int(bin_ms)}ms_sse"] = float(np.sum(diff * diff))
        metrics[f"psth_{int(bin_ms)}ms_rmse"] = float(np.sqrt(np.mean(diff * diff))) if diff.size else np.nan
        metrics[f"psth_{int(bin_ms)}ms_pearson_r"] = safe_pearson(model_psth, data_psth)

    model_bins_20, data_bins_20 = binned_cache[20]
    model_bins_10, data_bins_10 = binned_cache[10]
    metrics["rcorr"] = safe_pearson(model_bins_10.sum(axis=0), data_bins_10.sum(axis=0))
    metrics["psth_20ms_cosine_similarity"] = safe_cosine(model_bins_20.sum(axis=0), data_bins_20.sum(axis=0))
    metrics["split_half_data_20ms"] = split_half_corr(data_bins_20)
    metrics["split_half_model_20ms"] = split_half_corr(model_bins_20)
    metrics["noise_corrected_r_10ms"] = noise_corrected_r(*binned_cache[10])
    metrics["lifetime_sparseness_data_20ms"] = lifetime_sparseness(data_bins_20.mean(axis=0))
    metrics["lifetime_sparseness_model_20ms"] = lifetime_sparseness(model_bins_20.mean(axis=0))
    metrics["lifetime_sparseness_abs_error_20ms"] = abs(
        metrics["lifetime_sparseness_model_20ms"] - metrics["lifetime_sparseness_data_20ms"]
    )

    metrics["firing_rate_abs_error_hz"] = abs(mean_rate_hz(model_trial_time, dt_ms) - mean_rate_hz(data_trial_time, dt_ms))

    stimulus_info = get_stimulus_event_info(stimulus_path, response_bin_ms=20.0, n_bins=model_bins_20.shape[1])
    model_features = compute_response_features(model_bins_20, 20.0, stimulus_info)
    data_features = compute_response_features(data_bins_20, 20.0, stimulus_info)
    metrics["peak_rate_abs_error_hz"] = abs(model_features["peak_rate_hz"] - data_features["peak_rate_hz"])
    metrics["peak_latency_abs_error_ms"] = abs(model_features["peak_latency_ms"] - data_features["peak_latency_ms"])
    metrics["onset_responsiveness_abs_error"] = abs(
        model_features["onset_responsiveness"] - data_features["onset_responsiveness"]
    )
    metrics["offset_responsiveness_abs_error"] = abs(
        model_features["offset_responsiveness"] - data_features["offset_responsiveness"]
    )

    ordered_metrics = {
        "raw_raster_mse": metrics["raw_raster_mse"],
        "raw_raster_rmse": metrics["raw_raster_rmse"],
        "raw_raster_mae": metrics["raw_raster_mae"],
        "psth_5ms_sse": metrics["psth_5ms_sse"],
        "psth_10ms_sse": metrics["psth_10ms_sse"],
        "psth_20ms_sse": metrics["psth_20ms_sse"],
        "psth_50ms_sse": metrics["psth_50ms_sse"],
        "psth_5ms_pearson_r": metrics["psth_5ms_pearson_r"],
        "psth_10ms_pearson_r": metrics["psth_10ms_pearson_r"],
        "psth_20ms_pearson_r": metrics["psth_20ms_pearson_r"],
        "psth_50ms_pearson_r": metrics["psth_50ms_pearson_r"],
        "rcorr": metrics["rcorr"],
        "psth_20ms_cosine_similarity": metrics["psth_20ms_cosine_similarity"],
        "split_half_data_20ms": metrics["split_half_data_20ms"],
        "split_half_model_20ms": metrics["split_half_model_20ms"],
        "noise_corrected_r_10ms": metrics["noise_corrected_r_10ms"],
        "lifetime_sparseness_data_20ms": metrics["lifetime_sparseness_data_20ms"],
        "lifetime_sparseness_model_20ms": metrics["lifetime_sparseness_model_20ms"],
        "lifetime_sparseness_abs_error_20ms": metrics["lifetime_sparseness_abs_error_20ms"],
        "firing_rate_abs_error_hz": metrics["firing_rate_abs_error_hz"],
        "peak_rate_abs_error_hz": metrics["peak_rate_abs_error_hz"],
        "peak_latency_abs_error_ms": metrics["peak_latency_abs_error_ms"],
        "onset_responsiveness_abs_error": metrics["onset_responsiveness_abs_error"],
        "offset_responsiveness_abs_error": metrics["offset_responsiveness_abs_error"],
    }
    return ordered_metrics


def load_reference_joint_mds(joint_mds_dir: Path) -> Dict[str, np.ndarray]:
    summary_paths = [
        joint_mds_dir / "JointMDS_feature_summary.mat",
        joint_mds_dir / "SeparateMDS_feature_summary.mat",
    ]
    summary_path = next((path for path in summary_paths if path.exists()), None)
    if summary_path is None:
        raise FileNotFoundError(f"Could not find JointMDS summary MAT in {joint_mds_dir}")

    y_ref = np.asarray(loadmat(summary_path, variable_names=["Y"])["Y"], dtype=np.float64)

    output_bundle = loadmat(
        joint_mds_dir / "50epoch_100ms_and_300ms_loss_All_cells.mat",
        variable_names=["output", "losses"],
    )
    output = np.asarray(output_bundle["output"])
    losses = np.asarray(output_bundle["losses"])

    data_bundle = loadmat(
        joint_mds_dir / "all_units_info_with_polished_criteria_modified_perf.mat",
        variable_names=["all_data"],
        squeeze_me=True,
        struct_as_record=False,
    )
    all_data = data_bundle["all_data"]

    n_cells = int(output.shape[0])
    n_trials = int(output.shape[2])
    time_len = int(output.shape[-1])
    bin_width = 200

    final_loss = np.asarray(losses[-1, 1, :, :], dtype=np.float64)
    if final_loss.shape[0] != n_cells and final_loss.shape[1] == n_cells:
        final_loss = final_loss.T
    best_batch_idx = np.argmin(final_loss, axis=1)

    psth_data = np.zeros((n_cells, time_len // bin_width), dtype=np.float32)
    psth_model = np.zeros((n_cells, time_len // bin_width), dtype=np.float32)

    for cell_idx in range(n_cells):
        sim_raster = np.asarray(output[cell_idx, best_batch_idx[cell_idx]])
        sim_raster = squeeze_sim_raster(sim_raster, n_trials)
        psth_model[cell_idx] = summed_psth(sim_raster, bin_width)

        unit = all_data[cell_idx]
        data_raster = extract_unit_raster(unit, n_trials=n_trials, time_len=time_len)
        psth_data[cell_idx] = summed_psth(data_raster, bin_width)

    psths_all = np.vstack([psth_data, psth_model]).astype(np.float64)
    return {
        "y_ref": y_ref,
        "psths_all": psths_all,
        "n_cells": n_cells,
        "best_batch_idx": best_batch_idx.astype(int),
    }


def squeeze_sim_raster(sim_raster: np.ndarray, n_trials: int) -> np.ndarray:
    sim_raster = np.asarray(sim_raster, dtype=np.float32).squeeze()
    if sim_raster.ndim == 1:
        sim_raster = sim_raster.reshape(n_trials, -1)
    elif sim_raster.ndim == 3:
        if sim_raster.shape[0] == n_trials:
            sim_raster = sim_raster.sum(axis=1)
        else:
            sim_raster = sim_raster.reshape(n_trials, -1)
    elif sim_raster.ndim == 2 and sim_raster.shape[0] != n_trials and sim_raster.shape[1] == n_trials:
        sim_raster = sim_raster.T
    elif sim_raster.ndim != 2:
        sim_raster = sim_raster.reshape(n_trials, -1)
    return sim_raster


def extract_unit_raster(unit, *, n_trials: int, time_len: int) -> np.ndarray:
    focus = get_focus_index(getattr(unit, "tuning_type", "contra"))
    spike_times = unit.ctrl_tar1_timestamps
    if isinstance(spike_times, np.ndarray) and spike_times.ndim == 2:
        focus = min(max(focus, 0), spike_times.shape[1] - 1)
        spike_times = spike_times[:, focus]

    raster = np.zeros((n_trials, time_len), dtype=np.float32)
    for trial_idx in range(n_trials):
        trial_times = np.asarray(spike_times[trial_idx], dtype=np.float64).reshape(-1)
        if trial_times.size == 0:
            continue
        mask = np.isfinite(trial_times) & (trial_times > 0.0) & (trial_times < 2.9801)
        sample_idx = np.round(trial_times[mask] * 10000.0).astype(int)
        sample_idx = sample_idx[(sample_idx >= 0) & (sample_idx < time_len)]
        raster[trial_idx, sample_idx] = 1.0
    return raster


def get_focus_index(tuning_type) -> int:
    tuning_text = str(tuning_type).lower()
    if "contra" in tuning_text:
        return 0
    if "45" in tuning_text:
        return 1
    if "center" in tuning_text:
        return 2
    if "ipsi" in tuning_text:
        return 3
    return 0


def project_point_to_reference_space(candidate_psth: np.ndarray, reference_psths: np.ndarray, y_ref: np.ndarray) -> np.ndarray:
    candidate_psth = np.asarray(candidate_psth, dtype=np.float64)
    dists = np.linalg.norm(reference_psths - candidate_psth[None, :], axis=1)
    start_xy = np.asarray(y_ref[int(np.argmin(dists))], dtype=np.float64)

    def residuals(xy: np.ndarray) -> np.ndarray:
        return np.linalg.norm(y_ref - xy[None, :], axis=1) - dists

    result = least_squares(residuals, x0=start_xy, method="trf")
    return result.x.astype(np.float64)


def save_report_figure(output_path: Path, metrics: Dict[str, float], projection_payload: Dict[str, object]) -> None:
    y_ref = np.asarray(projection_payload["reference_y"], dtype=np.float64)
    n_cells = int(projection_payload["n_cells"])
    reference_data_xy = np.asarray(projection_payload["reference_data_xy"], dtype=np.float64)
    reference_model_xy = np.asarray(projection_payload["reference_model_xy"], dtype=np.float64)
    manual_fit_xy = np.asarray(projection_payload["manual_fit_xy"], dtype=np.float64)
    unit_index = int(projection_payload["unit_index"])
    best_batch_idx = int(projection_payload["best_batch_idx"])
    reference_best_batch_idx = int(projection_payload["reference_best_batch_idx"])

    limits = compute_square_limits(
        np.vstack([y_ref, reference_data_xy[None, :], reference_model_xy[None, :], manual_fit_xy[None, :]])
    )

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 1.0], height_ratios=[1.0, 1.0])
    ax_joint = fig.add_subplot(gs[0, 0])
    ax_data = fig.add_subplot(gs[0, 1])
    ax_model = fig.add_subplot(gs[1, 0])
    ax_text = fig.add_subplot(gs[1, 1])

    draw_reference_joint(ax_joint, y_ref, n_cells)
    ax_joint.plot(
        [reference_data_xy[0], reference_model_xy[0]],
        [reference_data_xy[1], reference_model_xy[1]],
        color="#808080",
        linewidth=1.0,
        alpha=0.8,
        zorder=4,
    )
    ax_joint.plot(
        [reference_data_xy[0], manual_fit_xy[0]],
        [reference_data_xy[1], manual_fit_xy[1]],
        color="#b3471f",
        linewidth=1.0,
        linestyle="--",
        alpha=0.9,
        zorder=4,
    )
    ax_joint.scatter(
        reference_data_xy[0],
        reference_data_xy[1],
        s=150,
        color="#1f77b4",
        edgecolor="black",
        linewidth=0.9,
        zorder=6,
        label=f"Cell {unit_index} data",
    )
    ax_joint.scatter(
        reference_model_xy[0],
        reference_model_xy[1],
        s=150,
        marker="s",
        color="#ff9f1c",
        edgecolor="black",
        linewidth=0.9,
        zorder=6,
        label=f"Cell {unit_index} original model",
    )
    ax_joint.scatter(
        manual_fit_xy[0],
        manual_fit_xy[1],
        s=220,
        marker="*",
        color="#d62728",
        edgecolor="black",
        linewidth=0.9,
        zorder=7,
        label=f"Cell {unit_index} manual fit",
    )
    annotate_projection_points(ax_joint, unit_index, reference_data_xy, reference_model_xy, manual_fit_xy)
    ax_joint.set_title(
        f"Joint MDS projection for unit {unit_index}\n"
        f"(reference batch = {reference_best_batch_idx}, manual best batch = {best_batch_idx})"
    )
    ax_joint.legend(loc="upper right", fontsize=8)
    apply_axes(ax_joint, limits)

    ax_data.scatter(y_ref[:n_cells, 0], y_ref[:n_cells, 1], s=18, color="#4f6db8", alpha=0.65)
    ax_data.scatter(reference_data_xy[0], reference_data_xy[1], s=150, color="#1f77b4", edgecolor="black", linewidth=0.9)
    ax_data.scatter(manual_fit_xy[0], manual_fit_xy[1], s=220, marker="*", color="#d62728", edgecolor="black", linewidth=0.9)
    annotate_projection_points(ax_data, unit_index, reference_data_xy, None, manual_fit_xy)
    ax_data.set_title("Reference data cloud + cell data/manual fit")
    apply_axes(ax_data, limits)

    ax_model.scatter(y_ref[n_cells:, 0], y_ref[n_cells:, 1], s=20, color="#c97f2f", alpha=0.65)
    ax_model.scatter(reference_model_xy[0], reference_model_xy[1], s=150, marker="s", color="#ff9f1c", edgecolor="black", linewidth=0.9)
    ax_model.scatter(manual_fit_xy[0], manual_fit_xy[1], s=220, marker="*", color="#d62728", edgecolor="black", linewidth=0.9)
    annotate_projection_points(ax_model, unit_index, None, reference_model_xy, manual_fit_xy)
    ax_model.set_title("Reference model cloud + original/manual fit")
    apply_axes(ax_model, limits)

    ax_text.axis("off")
    metric_lines = [f"{name:>32}: {value: .6g}" for name, value in metrics.items()]
    ax_text.text(
        0.0,
        1.0,
        "\n".join(metric_lines),
        va="top",
        ha="left",
        family="monospace",
        fontsize=9,
    )
    ax_text.set_title(f"{len(metrics)} post-fit metrics", loc="left")

    fig.tight_layout()
    fig.savefig(output_path, format="svg", bbox_inches="tight")
    plt.close(fig)


def draw_reference_joint(ax, y_ref: np.ndarray, n_cells: int) -> None:
    for cell_idx in range(n_cells):
        ax.plot(
            [y_ref[cell_idx, 0], y_ref[n_cells + cell_idx, 0]],
            [y_ref[cell_idx, 1], y_ref[n_cells + cell_idx, 1]],
            color="#d0d0d0",
            linewidth=0.5,
            zorder=1,
        )
    ax.scatter(y_ref[:n_cells, 0], y_ref[:n_cells, 1], s=18, color="#4f6db8", alpha=0.35)
    ax.scatter(y_ref[n_cells:, 0], y_ref[n_cells:, 1], s=20, color="#c97f2f", alpha=0.35, marker="s")


def annotate_projection_points(
    ax,
    unit_index: int,
    reference_data_xy: np.ndarray | None,
    reference_model_xy: np.ndarray | None,
    manual_fit_xy: np.ndarray | None,
) -> None:
    if reference_data_xy is not None:
        ax.annotate(
            f"Cell {unit_index} data",
            xy=reference_data_xy,
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=8,
            color="#1f77b4",
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.8),
        )
    if reference_model_xy is not None:
        ax.annotate(
            f"Cell {unit_index} original",
            xy=reference_model_xy,
            xytext=(6, -12),
            textcoords="offset points",
            fontsize=8,
            color="#9a5a00",
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.8),
        )
    if manual_fit_xy is not None:
        ax.annotate(
            f"Cell {unit_index} manual",
            xy=manual_fit_xy,
            xytext=(6, 10),
            textcoords="offset points",
            fontsize=8,
            color="#a61c1c",
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.8),
        )


def compute_square_limits(points: np.ndarray, padding_frac: float = 0.05) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    x = points[:, 0]
    y = points[:, 1]
    x_mid = 0.5 * (np.nanmin(x) + np.nanmax(x))
    y_mid = 0.5 * (np.nanmin(y) + np.nanmax(y))
    half_span = 0.5 * max(np.nanmax(x) - np.nanmin(x), np.nanmax(y) - np.nanmin(y))
    half_span *= 1.0 + padding_frac
    return (x_mid - half_span, x_mid + half_span), (y_mid - half_span, y_mid + half_span)


def apply_axes(ax, limits: Tuple[Tuple[float, float], Tuple[float, float]]) -> None:
    ax.set_xlim(limits[0])
    ax.set_ylim(limits[1])
    ax.set_xlabel("MDS 1")
    ax.set_ylabel("MDS 2")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, alpha=0.25)


def write_metric_csv(path: Path, metrics: Dict[str, float]) -> None:
    rows = [{"metric": key, "value": value} for key, value in metrics.items()]
    write_dict_rows_csv(path, rows)


def write_dict_rows_csv(path: Path, rows: Iterable[Dict[str, object]]) -> None:
    rows = list(rows)
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    headers: List[str] = list(rows[0].keys())
    lines = [",".join(headers)]
    for row in rows:
        values = []
        for header in headers:
            value = row.get(header, "")
            if isinstance(value, float):
                values.append(f"{value:.12g}")
            else:
                values.append(str(value))
        lines.append(",".join(values))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
