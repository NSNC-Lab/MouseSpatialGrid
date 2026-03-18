from __future__ import annotations

import argparse
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.io.matlab import mat_struct
from scipy.io.wavfile import read as wavread


def select_output_mat(default_dir: Path) -> Path:
    mat_files = sorted(default_dir.glob("output_compressed_Eprop*.mat"))
    if not mat_files:
        raise FileNotFoundError(f"No output_compressed_Eprop*.mat files in {default_dir}")

    try:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        root.update()
        selected = filedialog.askopenfilename(
            title="Select an output .mat file",
            initialdir=str(default_dir),
            filetypes=[("MAT files", "*.mat")],
        )
        root.destroy()
        if selected:
            return Path(selected)
    except Exception:
        pass

    return mat_files[-1]


def infer_unit_index_from_filename(mat_path: Path, fallback: int) -> int:
    match = re.search(r"Cell(\d+)", mat_path.stem)
    if match:
        return int(match.group(1))
    return fallback


def _to_trial_time(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr)

    if arr.ndim == 2:
        return arr

    if arr.ndim == 3:
        candidate_axes = [axis for axis, size in enumerate(arr.shape) if size == 10]
        if not candidate_axes:
            raise ValueError(f"Could not find trial axis of length 10 in shape {arr.shape}")

        trial_axis = candidate_axes[0]
        moved = np.moveaxis(arr, trial_axis, 0)
        moved = np.squeeze(moved)

        if moved.ndim == 1:
            moved = moved[None, :]
        if moved.ndim != 2:
            raise ValueError(f"Cannot coerce array with shape {arr.shape} to (trials,time)")
        return moved

    if arr.ndim == 4:
        # Expected: (batch, trials, channels, time) or (batch, trials, time, channels)
        best = arr[0]
        best = np.squeeze(best)
        if best.ndim == 2:
            return best

        if best.ndim == 3:
            candidate_axes = [axis for axis, size in enumerate(best.shape) if size == 10]
            if not candidate_axes:
                raise ValueError(f"Cannot infer trial axis from shape {best.shape}")
            trial_axis = candidate_axes[0]
            moved = np.moveaxis(best, trial_axis, 0)
            moved = np.squeeze(moved)
            if moved.ndim == 2:
                return moved

    raise ValueError(f"Unsupported output shape: {arr.shape}")


def load_model_raster(output_mat: Path) -> np.ndarray:
    mat = loadmat(output_mat)

    if "best_output" in mat:
        model = mat["best_output"]
    elif "output" in mat:
        model = mat["output"]
    else:
        keys = [k for k in mat.keys() if not k.startswith("__")]
        raise KeyError(f"Could not find 'best_output' or 'output' in {output_mat}. Keys: {keys}")

    raster = _to_trial_time(model).astype(np.float64)

    if raster.shape[0] != 10:
        raise ValueError(f"Expected 10 trials, got shape {raster.shape}")

    raster = (raster > 0).astype(np.float64)
    return raster


def load_target1_raster(all_units_path: Path, unit_index_1based: int, timesteps: int, dt_ms: float) -> np.ndarray:
    mat = loadmat(all_units_path, variable_names=["all_data"], squeeze_me=True, struct_as_record=False)
    all_data = mat["all_data"]

    if not (1 <= unit_index_1based <= len(all_data)):
        raise IndexError(
            f"unit_index out of range: {unit_index_1based}, available 1..{len(all_data)}"
        )

    unit = all_data[unit_index_1based - 1]
    if not isinstance(unit, mat_struct):
        raise TypeError("Unexpected MATLAB struct format for selected unit")

    spike_times = unit.ctrl_tar1_timestamps
    spike_times = np.asarray(spike_times)
    if spike_times.ndim == 2:
        spike_times = spike_times[:, 0]

    target_raster = np.zeros((10, timesteps), dtype=np.float64)
    duration_s = timesteps * dt_ms / 1000.0

    for trial in range(10):
        times = np.asarray(spike_times[trial]).squeeze()
        times = times[(times >= 0) & (times < duration_s)]
        inds = np.rint(times * (1000.0 / dt_ms)).astype(int)
        inds = inds[(inds >= 0) & (inds < timesteps)]
        target_raster[trial, inds] = 1.0

    return target_raster


def psth_from_raster(raster: np.ndarray, bin_width_ms: float, dt_ms: float) -> np.ndarray:
    bin_steps = int(round(bin_width_ms / dt_ms))
    if bin_steps < 1:
        raise ValueError(f"bin_width_ms={bin_width_ms} is too small for dt_ms={dt_ms}")

    timesteps = raster.shape[1]
    usable = (timesteps // bin_steps) * bin_steps
    if usable == 0:
        raise ValueError("Not enough timesteps for one PSTH bin")

    reshaped = raster[:, :usable].reshape(raster.shape[0], -1, bin_steps)
    counts = reshaped.sum(axis=2)
    return counts.mean(axis=0)


def waveform_for_target1(path: Path, stim_key: str = "target1") -> np.ndarray:
    if path.suffix.lower() == ".mat":
        stim_mat = loadmat(path, squeeze_me=True, struct_as_record=False)
        if stim_key not in stim_mat:
            keys = [k for k in stim_mat.keys() if not k.startswith("__")]
            raise KeyError(f"Stim key '{stim_key}' not found in {path}. Available keys: {keys}")
        audio = np.asarray(stim_mat[stim_key], dtype=np.float64)
    else:
        _, audio = wavread(path)
        audio = np.asarray(audio, dtype=np.float64)

    if audio.ndim == 2:
        audio = audio.mean(axis=1)

    peak = np.max(np.abs(audio))
    if peak > 0:
        audio = audio / peak

    return audio


def plot_layout(
    target_raster: np.ndarray,
    model_raster: np.ndarray,
    stim_waveform: np.ndarray,
    dt_ms: float,
    psth_bin_ms: float,
    unit_index: int,
    mat_name: str,
    font_size: int = 16,
) -> None:
    target_psth = psth_from_raster(target_raster, psth_bin_ms, dt_ms)
    model_psth = psth_from_raster(model_raster, psth_bin_ms, dt_ms)
    duration_s = model_raster.shape[1] * dt_ms / 1000.0

    # Compute cross-correlation
    if len(target_psth) == len(model_psth) and len(target_psth) > 0:
        # Normalize the signals
        target_norm = (target_psth - np.mean(target_psth)) / (np.std(target_psth) + 1e-10)
        model_norm = (model_psth - np.mean(model_psth)) / (np.std(model_psth) + 1e-10)
        
        # Compute normalized cross-correlation
        xcorr = np.correlate(target_norm, model_norm, mode='full') / len(target_psth)
        max_xcorr = np.max(xcorr)
        lag_at_max = np.argmax(xcorr) - (len(target_psth) - 1)
        
        print(f"PSTH Cross-correlation: max={max_xcorr:.4f} at lag={lag_at_max} bins")
    else:
        print(f"Warning: PSTH length mismatch - target: {len(target_psth)}, model: {len(model_psth)}")

    # Configure global font sizes
    plt.rcParams.update({
        'font.size': font_size,
        'axes.labelsize': font_size,
        'axes.titlesize': font_size,
        'xtick.labelsize': font_size,
        'ytick.labelsize': font_size,
        'legend.fontsize': font_size,
        'figure.titlesize': font_size + 2,
    })

    fig, axes = plt.subplots(4, 1, figsize=(10, 8), gridspec_kw={"height_ratios": [1.2, 1.0, 1.0, 1.2]})

    # 1) PSTH (top)
    ax = axes[0]
    x = np.arange(target_psth.shape[0])
    ax.plot(x, target_psth, 'k', label="target1 PSTH", linewidth=1.2)
    ax.plot(x, model_psth, label="model PSTH", linewidth=1.2)
    #ax.set_ylabel("Spikes/bin")
    ax.set_xticks([])  
    ax.set_yticks([])
    ax.set_xlim([0, max(1, target_psth.shape[0] - 1)])
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.25)

    # 2) Target raster
    ax = axes[1]
    target_events = [np.flatnonzero(target_raster[tr] > 0) for tr in range(target_raster.shape[0])]
    ax.eventplot(target_events, linewidths=0.8,color='black')
    ax.set_ylabel("Target")
    #ax.set_title("Target 1 raster")
    ax.set_xticks([])  
    ax.set_yticks([])
    ax.set_xlim([0, target_raster.shape[1]])
    ax.set_ylim([0, target_raster.shape[0]])

    # 3) Model raster
    ax = axes[2]
    model_events = [np.flatnonzero(model_raster[tr] > 0) for tr in range(model_raster.shape[0])]
    ax.eventplot(model_events, linewidths=0.8)
    ax.set_ylabel("Model")
    #ax.set_title("Model raster")
    ax.set_xticks([])  
    ax.set_yticks([])
    ax.set_xlim([0, model_raster.shape[1]])
    ax.set_ylim([0, model_raster.shape[0]])

    # 4) Stimulus waveform (bottom)

    #stim_waveform = stim_waveform[62000:]

    ax = axes[3]
    stim_time_s = np.linspace(0.0, duration_s, stim_waveform.shape[0], endpoint=False)
    ax.plot(stim_time_s, stim_waveform, linewidth=0.8)
    #ax.set_ylabel("Norm amp")
    ax.set_xlabel("Time (s)")
    #ax.set_title("Target 1 stimulus")
    ax.set_yticks([])
    ax.grid(True, alpha=0.2)
    ax.set_xlim([0.0, duration_s])
    stim_ticks = np.linspace(0.0, duration_s, 7)
    ax.set_xticks(stim_ticks)
    ax.set_xticklabels([f"{tick:.2f}" for tick in stim_ticks])
    

    fig.suptitle(f"Unit {unit_index}")
    plt.tight_layout()
    save_name = f"psth_raster_stim_unit{unit_index}_{Path(mat_name).stem}.svg"
    fig.savefig(save_name, format="svg", bbox_inches="tight")
    print(f"Saved SVG: {save_name}")
    plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot PSTH + rasters + stim for selected E-prop output")
    parser.add_argument(
        "--mat-file",
        type=str,
        default="",
        help="Optional explicit output .mat path; if omitted, opens picker in E-prop directory.",
    )
    parser.add_argument(
        "--eprop-dir",
        type=str,
        default=r"C:\Users\ipboy\Documents\GitHub\MouseSpatialGrid\LearningModels\E-prop",
        help="Directory containing output_compressed_Eprop*.mat files.",
    )
    parser.add_argument(
        "--all-units-mat",
        type=str,
        default=r"C:\Users\ipboy\Documents\GitHub\ModelingEffort\Multi-Channel\Plotting\OliverDataPlotting\all_units_info_with_polished_criteria_modified_perf.mat",
        help="Path to all_units_info_with_polished_criteria_modified_perf.mat",
    )
    parser.add_argument(
        "--unit-index",
        type=int,
        default=200,
        help="1-based unit index for all_data (defaults to 200 to match run_main_integrated.py).",
    )
    parser.add_argument(
        "--dt-ms",
        type=float,
        default=0.1,
        help="Simulation dt in ms (used for target spike binning and PSTH binning).",
    )
    parser.add_argument(
        "--psth-bin-ms",
        type=float,
        default=20.0,
        help="PSTH bin width in ms.",
    )
    parser.add_argument(
        "--stim-wav",
        type=str,
        default=r"C:\Users\ipboy\Documents\GitHub\ModelingEffort\Multi-Channel\Plotting\OliverDataPlotting\sound_files.mat",
        help="Stimulus source for bottom panel (.mat or .wav).",
    )
    parser.add_argument(
        "--stim-key",
        type=str,
        default="target1",
        help="Variable name inside --stim-wav when it is a .mat file.",
    )
    parser.add_argument(
        "--font-size",
        type=int,
        default=16,
        help="Base font size for all text in the plot (default: 16).",
    )

    args = parser.parse_args()

    eprop_dir = Path(args.eprop_dir)
    if args.mat_file:
        output_mat = Path(args.mat_file)
    else:
        output_mat = select_output_mat(eprop_dir)

    if not output_mat.exists():
        raise FileNotFoundError(f"Selected output mat does not exist: {output_mat}")

    unit_index = infer_unit_index_from_filename(output_mat, args.unit_index)

    model_raster = load_model_raster(output_mat)
    target_raster = load_target1_raster(
        Path(args.all_units_mat),
        unit_index,
        timesteps=model_raster.shape[1],
        dt_ms=args.dt_ms,
    )
    stim_waveform = waveform_for_target1(Path(args.stim_wav), stim_key=args.stim_key)

    plot_layout(
        target_raster=target_raster,
        model_raster=model_raster,
        stim_waveform=stim_waveform,
        dt_ms=args.dt_ms,
        psth_bin_ms=args.psth_bin_ms,
        unit_index=unit_index,
        mat_name=output_mat.name,
        font_size=args.font_size,
    )


if __name__ == "__main__":
    main()
