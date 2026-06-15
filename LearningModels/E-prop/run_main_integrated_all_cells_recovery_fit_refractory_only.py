import set_options
import declarations
import declarations_output_inhibit
from BuildFile import Forwards_Method_cupy_recovery_batch
import numpy as np
import time
from datetime import datetime
from scipy.io import loadmat, savemat
import yaml
import os
from pathlib import Path
import importlib.util
import sys
from strf_handler import call_strfs
from input_handler import call_inputs
import matplotlib.pyplot as plt


CHECKPOINT_EVERY_EPOCHS = 1
CHECKPOINT_PREFIX = "checkpoint_Eprop_All_cells_recovery_fit_refractory_only"
TOTAL_EPOCHS = 45
RESUME_CHECKPOINT_PATH = ""
START_PARAMS_PATH = Path(
    r"C:\Users\ipboy\Documents\GitHub\ModelingEffort\Multi-Channel\Plotting\OliverDataPlotting\checkpoint_Eprop_All_cells_recovery_fit_refractory_only_epoch_0020.mat"
)
RECOVERY_FUNCS_PATH = Path(r"C:\Users\ipboy\Desktop\SingleChannelFigures\Figure1_MDS_Plots\recovery_funcs_probabilistic.mat")
RECOVERY_FUNCS_KEY = "ws"
GENERATED_RECOVERY_SOLVER_PATH = Path(__file__).resolve().parent / "BuildFile" / "generated_solve_file_recovery_fit_refractory_only.py"
FIX_INPUTS_ACROSS_REFRACTORY_CANDIDATES = True

# The refractory gate is stochastic/discrete, so this runner learns it with a
# batch-wise evolutionary search rather than backpropagating through the random
# draw. The learned recovery curve is:
#
#   p_allow(t) = ceiling * (1 - exp(-((t - abs_ref) / tau) ** shape))
#
# with p_allow(t <= abs_ref) = 0. This is a hazard-style per-threshold-crossing
# gate, not a CDF of observed ISIs.
REFRACTORY_RANDOM_SEED = 17
REFRACTORY_MIN_CEILING = 0.05
REFRACTORY_PARAM_BOUNDS = {
    "log_abs_steps": (np.log(1.0), np.log(100.0)),
    "log_tau_steps": (np.log(2.0), np.log(5000.0)),
    "log_shape": (np.log(0.7), np.log(8.0)),
    "ceiling_raw": (-4.0, 6.0),
}
REFRACTORY_INITIAL_SIGMA = np.asarray([0.65, 1.00, 0.45, 1.50], dtype=np.float32)
REFRACTORY_MIN_SIGMA = np.asarray([0.08, 0.10, 0.06, 0.20], dtype=np.float32)
REFRACTORY_SIGMA_DECAY = 0.88
REFRACTORY_CENTER_MOMENTUM = 0.35

REFRACTORY_ISI_HIST_WEIGHT = 0.5
REFRACTORY_ISI_MEDIAN_WEIGHT = 1
REFRACTORY_ISI_COUNT_WEIGHT = 0.25
REFRACTORY_ISI_CV_WEIGHT = 1


def load_recovery_functions(path=RECOVERY_FUNCS_PATH, key=RECOVERY_FUNCS_KEY, num_cells=None):
    recovery_path = Path(path)
    if not recovery_path.exists():
        raise FileNotFoundError(f"Recovery function file not found: {recovery_path}")

    mat = loadmat(str(recovery_path), variable_names=[key], squeeze_me=True, struct_as_record=False)
    if key not in mat:
        available = [k for k in mat.keys() if not k.startswith("__")]
        raise KeyError(f"Recovery function key '{key}' not found in {recovery_path}. Available keys: {available}")

    recovery_funcs = np.asarray(mat[key], dtype=np.float32)
    if recovery_funcs.ndim == 1:
        recovery_funcs = recovery_funcs.reshape(1, -1)

    if num_cells is not None and recovery_funcs.shape[0] != num_cells and recovery_funcs.shape[1] == num_cells:
        recovery_funcs = recovery_funcs.T

    if num_cells is not None and recovery_funcs.shape[0] != num_cells:
        raise ValueError(
            f"Recovery functions should have one row per cell. "
            f"Expected {num_cells} rows, got shape {recovery_funcs.shape}."
        )

    recovery_funcs = np.nan_to_num(recovery_funcs, nan=1.0, posinf=1.0, neginf=0.0)
    return np.clip(recovery_funcs, 0.0, 1.0).astype(np.float32, copy=False)


def _coerce_per_cell_params(params, num_params, num_cells):
    params = np.asarray(params, dtype=np.float64)
    if params.shape == (num_params, num_cells):
        return params
    if params.shape == (num_cells, num_params):
        return params.T
    raise ValueError(
        f"Expected per-cell params with shape {(num_params, num_cells)} "
        f"or {(num_cells, num_params)}, got {params.shape}."
    )


def load_fixed_params_from_mat(path, batch_size, num_params, num_cells):
    params_path = Path(path)
    if not params_path.exists():
        raise FileNotFoundError(f"Starting parameter file not found: {params_path}")

    checkpoint = loadmat(str(params_path), squeeze_me=True, struct_as_record=False)
    if "best_params_per_cell" in checkpoint:
        fixed_params_per_cell = _coerce_per_cell_params(
            checkpoint["best_params_per_cell"],
            num_params,
            num_cells,
        )
    elif "params" in checkpoint:
        params = np.asarray(checkpoint["params"], dtype=np.float64)
        if params.shape != (num_params, num_cells, batch_size):
            raise ValueError(
                f"Expected params with shape {(num_params, num_cells, batch_size)}, "
                f"got {params.shape}."
            )

        if "best_batch_id_per_cell" in checkpoint:
            best_batch_id = np.asarray(checkpoint["best_batch_id_per_cell"]).reshape(num_cells).astype(np.int64) - 1
            best_batch_id = np.clip(best_batch_id, 0, batch_size - 1)
            fixed_params_per_cell = params[:, np.arange(num_cells), best_batch_id]
        else:
            fixed_params_per_cell = params[:, :, 0]
    else:
        available = [key for key in checkpoint.keys() if not key.startswith("__")]
        raise KeyError(f"No params or best_params_per_cell found in {params_path}. Available keys: {available}")

    p = np.repeat(fixed_params_per_cell[:, :, None], batch_size, axis=2)
    return p.astype(np.float64, copy=False), fixed_params_per_cell.astype(np.float32, copy=False)


def load_start_refractory_params_from_mat(path, num_cells):
    params_path = Path(path)
    if not params_path.exists():
        return None

    checkpoint = loadmat(
        str(params_path),
        variable_names=["best_refractory_params_per_cell"],
        squeeze_me=True,
        struct_as_record=False,
    )
    if "best_refractory_params_per_cell" not in checkpoint:
        return None

    refractory_params = np.asarray(checkpoint["best_refractory_params_per_cell"], dtype=np.float32)
    if refractory_params.shape == (4, num_cells):
        refractory_params = refractory_params.T
    if refractory_params.shape != (num_cells, 4):
        raise ValueError(
            f"Expected best_refractory_params_per_cell with shape {(num_cells, 4)}, "
            f"got {refractory_params.shape}."
        )
    return clip_refractory_params(refractory_params)


def repeat_first_batch(array, batch_axis, batch_size):
    array = np.asarray(array)
    if array.ndim <= batch_axis or array.shape[batch_axis] != batch_size:
        return array
    first_batch = np.take(array, indices=[0], axis=batch_axis)
    return np.repeat(first_batch, batch_size, axis=batch_axis)


def logit(x):
    x = np.clip(x, 1e-6, 1.0 - 1e-6)
    return np.log(x / (1.0 - x))


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def clip_refractory_params(params):
    params = np.asarray(params, dtype=np.float32).copy()
    bounds = [
        REFRACTORY_PARAM_BOUNDS["log_abs_steps"],
        REFRACTORY_PARAM_BOUNDS["log_tau_steps"],
        REFRACTORY_PARAM_BOUNDS["log_shape"],
        REFRACTORY_PARAM_BOUNDS["ceiling_raw"],
    ]

    for param_idx, (lo, hi) in enumerate(bounds):
        params[..., param_idx] = np.clip(params[..., param_idx], lo, hi)

    return params


def refractory_params_to_values(params):
    params = clip_refractory_params(params)
    abs_steps = np.exp(params[..., 0])
    tau_steps = np.exp(params[..., 1])
    shape = np.exp(params[..., 2])
    ceiling = REFRACTORY_MIN_CEILING + (1.0 - REFRACTORY_MIN_CEILING) * sigmoid(params[..., 3])
    return abs_steps.astype(np.float32), tau_steps.astype(np.float32), shape.astype(np.float32), ceiling.astype(np.float32)


def estimate_initial_refractory_params(data_rasters):
    num_cells = data_rasters.shape[0]
    params = np.zeros((num_cells, 4), dtype=np.float32)

    for cell_idx in range(num_cells):
        isis = collect_isis_steps(data_rasters[cell_idx])
        if isis.size:
            abs_steps = np.percentile(isis, 1)
            tau_steps = max(np.percentile(isis, 15) - abs_steps, 2.0)
        else:
            abs_steps = 5.0
            tau_steps = 25.0

        abs_steps = np.clip(abs_steps, 1.0, 100.0)
        tau_steps = np.clip(tau_steps, 2.0, 5000.0)
        shape = 2.5
        ceiling = 0.98
        ceiling_norm = (ceiling - REFRACTORY_MIN_CEILING) / (1.0 - REFRACTORY_MIN_CEILING)

        params[cell_idx] = [
            np.log(abs_steps),
            np.log(tau_steps),
            np.log(shape),
            logit(ceiling_norm),
        ]

    return clip_refractory_params(params)


def propose_refractory_params(center_params, sigma, batch_size, rng):
    center_params = clip_refractory_params(center_params)
    sigma = np.maximum(np.asarray(sigma, dtype=np.float32), REFRACTORY_MIN_SIGMA)

    candidates = center_params[:, None, :] + rng.normal(
        loc=0.0,
        scale=sigma[None, None, :],
        size=(center_params.shape[0], batch_size, center_params.shape[1]),
    ).astype(np.float32)

    # Preserve one incumbent candidate per cell, so an epoch cannot lose the
    # current best refractory setting purely because of sampling noise.
    candidates[:, 0, :] = center_params
    return clip_refractory_params(candidates)


def build_parametric_recovery_functions(refractory_param_candidates, num_steps):
    abs_steps, tau_steps, shape, ceiling = refractory_params_to_values(refractory_param_candidates)
    x = np.arange(num_steps, dtype=np.float32)[None, None, :]

    elapsed = np.maximum(x - abs_steps[:, :, None], 0.0)
    tau = np.maximum(tau_steps[:, :, None], np.finfo(np.float32).eps)
    power = np.clip(shape[:, :, None], 0.2, 12.0)
    recovery = ceiling[:, :, None] * (1.0 - np.exp(-((elapsed / tau) ** power)))
    recovery[elapsed <= 0.0] = 0.0
    recovery = np.nan_to_num(recovery, nan=0.0, posinf=1.0, neginf=0.0)

    return np.clip(recovery, 0.0, 1.0).astype(np.float32, copy=False)


def collect_isis_steps(raster):
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
        raise ValueError(f"Expected 1D, 2D, or 3D raster, got shape {raster.shape}")

    if not isis:
        return np.empty((0,), dtype=np.float32)

    out = np.concatenate(isis).astype(np.float32, copy=False)
    return out[np.isfinite(out) & (out > 0)]


def make_isi_score_bin_edges(time_len):
    early_edges = np.arange(0, 501, 5, dtype=np.int32)
    late_edges = np.arange(550, time_len + 50, 50, dtype=np.int32)
    return np.unique(np.concatenate([early_edges, late_edges]))


def normalized_histogram(values, bin_edges):
    counts, _ = np.histogram(values, bins=bin_edges)
    total = counts.sum()
    if total <= 0:
        return np.zeros(len(bin_edges) - 1, dtype=np.float32)
    return (counts / total).astype(np.float32)


def precompute_data_isi_reference(data_rasters, bin_edges):
    num_cells = data_rasters.shape[0]
    hist = np.zeros((num_cells, len(bin_edges) - 1), dtype=np.float32)
    median_steps = np.full((num_cells,), np.nan, dtype=np.float32)
    cv = np.full((num_cells,), np.nan, dtype=np.float32)
    n_isis = np.zeros((num_cells,), dtype=np.int32)

    for cell_idx in range(num_cells):
        isis = collect_isis_steps(data_rasters[cell_idx])
        hist[cell_idx] = normalized_histogram(isis, bin_edges)
        n_isis[cell_idx] = isis.size
        if isis.size:
            median_steps[cell_idx] = np.median(isis)
        if isis.size >= 2 and np.mean(isis) > 0:
            cv[cell_idx] = np.std(isis) / np.mean(isis)

    return {"hist": hist, "median_steps": median_steps, "cv": cv, "n_isis": n_isis}


def safe_log_ratio(a, b):
    if not np.isfinite(a) or not np.isfinite(b):
        return 0.0
    return abs(np.log((float(a) + 1.0) / (float(b) + 1.0)))


def score_refractory_candidates(output, data_isi_reference, isi_bin_edges):
    output = np.asarray(output)
    num_cells, batch_size = output.shape[:2]

    scores = np.full((num_cells, batch_size), np.inf, dtype=np.float32)
    isi_hist_l1 = np.full((num_cells, batch_size), np.nan, dtype=np.float32)
    isi_median_penalty = np.full((num_cells, batch_size), np.nan, dtype=np.float32)
    isi_count_penalty = np.full((num_cells, batch_size), np.nan, dtype=np.float32)
    isi_cv_penalty = np.full((num_cells, batch_size), np.nan, dtype=np.float32)

    for cell_idx in range(num_cells):
        ref_hist = data_isi_reference["hist"][cell_idx]
        ref_median = data_isi_reference["median_steps"][cell_idx]
        ref_cv = data_isi_reference["cv"][cell_idx]
        ref_n = data_isi_reference["n_isis"][cell_idx]

        for batch_idx in range(batch_size):
            model_isis = collect_isis_steps(output[cell_idx, batch_idx])
            model_hist = normalized_histogram(model_isis, isi_bin_edges)

            hist_l1 = float(np.sum(np.abs(model_hist - ref_hist)))
            median_penalty = safe_log_ratio(np.median(model_isis) if model_isis.size else np.nan, ref_median)
            count_penalty = safe_log_ratio(model_isis.size, ref_n)
            model_cv = np.std(model_isis) / np.mean(model_isis) if model_isis.size >= 2 and np.mean(model_isis) > 0 else np.nan
            cv_penalty = safe_log_ratio(model_cv, ref_cv)

            isi_hist_l1[cell_idx, batch_idx] = hist_l1
            isi_median_penalty[cell_idx, batch_idx] = median_penalty
            isi_count_penalty[cell_idx, batch_idx] = count_penalty
            isi_cv_penalty[cell_idx, batch_idx] = cv_penalty
            scores[cell_idx, batch_idx] = (
                REFRACTORY_ISI_HIST_WEIGHT * hist_l1
                + REFRACTORY_ISI_MEDIAN_WEIGHT * median_penalty
                + REFRACTORY_ISI_COUNT_WEIGHT * count_penalty
                + REFRACTORY_ISI_CV_WEIGHT * cv_penalty
            )

    return scores, {
        "isi_hist_l1": isi_hist_l1,
        "isi_median_penalty": isi_median_penalty,
        "isi_count_penalty": isi_count_penalty,
        "isi_cv_penalty": isi_cv_penalty,
    }


def build_recovery_solver(solve_file_body, out_path=GENERATED_RECOVERY_SOLVER_PATH):
    solve_source = (
        "# pythran export solve_run(float64[:,:,:], float64[:,:,:], float64[:,:,:], float64[:,:,:], float64[:,:]) -> Tuple[float64[:,:,:,:], float64[:,:,:]]\n"
        "import numpy as np\n"
        "import cupy as cp\n"
        "from BuildFile import calculate_loss_eprop\n"
        "def solve_run(on_input,off_input,noise_token,rate_on,rate_off,rate_on_deriv,rate_off_deriv,data,p,recovery_funcs):\n"
        + solve_file_body
    )
    out_path = Path(out_path)
    out_path.write_text(solve_source, encoding="utf-8")
    print(f"{out_path.name} has been created at {out_path}.")
    module_name = "generated_solve_file_recovery_fit_refractory_only"
    spec = importlib.util.spec_from_file_location(module_name, out_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load generated recovery solver from {out_path}")
    solver = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = solver
    spec.loader.exec_module(solver)
    return solver

def save_training_checkpoint(
    epoch_one_indexed,
    p,
    m,
    v,
    t,
    losses,
    best_output_per_cell=None,
    best_loss_per_cell=None,
    best_batch_id_per_cell=None,
    best_epoch_per_cell=None,
    best_params_per_cell=None,
    best_refractory_params_per_cell=None,
    best_refractory_abs_steps_per_cell=None,
    best_refractory_tau_steps_per_cell=None,
    best_refractory_shape_per_cell=None,
    best_refractory_ceiling_per_cell=None,
    best_recovery_funcs_per_cell=None,
    refractory_param_history=None,
    refractory_score_history=None,
    best_isi_hist_l1_per_cell=None,
    best_isi_median_penalty_per_cell=None,
    best_isi_count_penalty_per_cell=None,
    best_isi_cv_penalty_per_cell=None,
    fixed_params_per_cell=None,
    prefix=CHECKPOINT_PREFIX,
):
    checkpoint_payload = {
        "epoch": np.asarray(epoch_one_indexed, dtype=np.int32),
        "params": np.asarray(p, dtype=np.float32),
        "adam_m": np.asarray(m, dtype=np.float32),
        "adam_v": np.asarray(v, dtype=np.float32),
        "adam_t": np.asarray(t, dtype=np.int32),
        "losses_so_far": np.asarray(losses, dtype=np.float32),
    }

    if best_output_per_cell is not None:
        checkpoint_payload.update(
            {
                "best_output_per_cell": np.asarray(best_output_per_cell, dtype=np.int8),
                "best_loss_per_cell": np.asarray(best_loss_per_cell, dtype=np.float32),
                "best_batch_id_per_cell": np.asarray(best_batch_id_per_cell, dtype=np.int32),
                "best_epoch_per_cell": np.asarray(best_epoch_per_cell, dtype=np.int32),
                "best_params_per_cell": np.asarray(best_params_per_cell, dtype=np.float32),
            }
        )

    optional_payload = {
        "best_refractory_params_per_cell": best_refractory_params_per_cell,
        "best_refractory_abs_steps_per_cell": best_refractory_abs_steps_per_cell,
        "best_refractory_tau_steps_per_cell": best_refractory_tau_steps_per_cell,
        "best_refractory_shape_per_cell": best_refractory_shape_per_cell,
        "best_refractory_ceiling_per_cell": best_refractory_ceiling_per_cell,
        "best_recovery_funcs_per_cell": best_recovery_funcs_per_cell,
        "refractory_param_history": refractory_param_history,
        "refractory_score_history": refractory_score_history,
        "best_isi_hist_l1_per_cell": best_isi_hist_l1_per_cell,
        "best_isi_median_penalty_per_cell": best_isi_median_penalty_per_cell,
        "best_isi_count_penalty_per_cell": best_isi_count_penalty_per_cell,
        "best_isi_cv_penalty_per_cell": best_isi_cv_penalty_per_cell,
        "fixed_params_per_cell": fixed_params_per_cell,
    }
    for key, value in optional_payload.items():
        if value is not None:
            checkpoint_payload[key] = np.asarray(value, dtype=np.float32)

    checkpoint_path = Path(f"{prefix}_epoch_{epoch_one_indexed:04d}.mat")
    latest_checkpoint_path = Path(f"{prefix}_latest.mat")
    latest_epoch_txt_path = Path(f"{prefix}_latest_epoch.txt")

    savemat(str(checkpoint_path), checkpoint_payload, do_compression=True)
    savemat(str(latest_checkpoint_path), checkpoint_payload, do_compression=True)
    latest_epoch_txt_path.write_text(f"{epoch_one_indexed}\n", encoding="utf-8")

    print(
        f"Checkpoint saved at epoch {epoch_one_indexed}: "
        f"{checkpoint_path} (latest: {latest_checkpoint_path})"
    )


def load_training_checkpoint(checkpoint_path, expected_shape):
    checkpoint_path = Path(checkpoint_path)
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")

    checkpoint = loadmat(str(checkpoint_path), squeeze_me=True, struct_as_record=False)

    required = ["epoch", "params", "adam_m", "adam_v", "adam_t"]
    missing = [key for key in required if key not in checkpoint]
    if missing:
        raise KeyError(
            f"Checkpoint missing required key(s): {missing}. "
            f"Available keys: {[k for k in checkpoint.keys() if not k.startswith('__')]}"
        )

    p_loaded = np.asarray(checkpoint["params"], dtype=np.float64)
    m_loaded = np.asarray(checkpoint["adam_m"], dtype=np.float64)
    v_loaded = np.asarray(checkpoint["adam_v"], dtype=np.float64)
    t_loaded = int(np.asarray(checkpoint["adam_t"]).squeeze())
    epoch_loaded = int(np.asarray(checkpoint["epoch"]).squeeze())

    if p_loaded.shape != expected_shape:
        raise ValueError(
            f"Checkpoint params shape {p_loaded.shape} does not match expected {expected_shape}."
        )
    if m_loaded.shape != expected_shape:
        raise ValueError(
            f"Checkpoint adam_m shape {m_loaded.shape} does not match expected {expected_shape}."
        )
    if v_loaded.shape != expected_shape:
        raise ValueError(
            f"Checkpoint adam_v shape {v_loaded.shape} does not match expected {expected_shape}."
        )

    best_state = {}
    optional_best_keys = [
        "best_output_per_cell",
        "best_loss_per_cell",
        "best_batch_id_per_cell",
        "best_epoch_per_cell",
        "best_params_per_cell",
        "best_refractory_params_per_cell",
        "best_refractory_abs_steps_per_cell",
        "best_refractory_tau_steps_per_cell",
        "best_refractory_shape_per_cell",
        "best_refractory_ceiling_per_cell",
        "best_recovery_funcs_per_cell",
        "best_isi_hist_l1_per_cell",
        "best_isi_median_penalty_per_cell",
        "best_isi_count_penalty_per_cell",
        "best_isi_cv_penalty_per_cell",
        "fixed_params_per_cell",
    ]
    for key in optional_best_keys:
        if key in checkpoint:
            best_state[key] = checkpoint[key]

    return p_loaded, m_loaded, v_loaded, t_loaded, epoch_loaded, best_state


class runSimulation(object):

    #Remove the outer loop later. Just for testing purposes.
    #for holder_thing in range(10):

    #Little control pannel for now (Eventually move this to a yaml file or whatever makes the most sense)

    gen_strfs_toggle = 0  #Toggle generating the STRFs
    gradients_toggle = 0  #Toggle generating the graidnets in the forwards process *Also toggles running epochs

    

    #Run STRF
    if gen_strfs_toggle == 1:
        call_strfs()
    #PreProcessesing  #Note! Recheck everything once you start running multichannel inputs -- check where the gain control for the tuning curves is and make sure we arn't doing extra steps
    #Will also need two worry about how exactly you are going to parse spks once you have multiple data streams

    num_cells = 220

    #num_cells, batch, trials, channels, timecourse
    #batch,trials,channels,timecourse

    #Set options
    opts = set_options.options()
    recovery_funcs = load_recovery_functions(num_cells=num_cells)
    print(f"Loaded recovery functions from {RECOVERY_FUNCS_PATH} with shape {recovery_funcs.shape}")
    #Declare architecture
    #arch = declarations.Declare_Architecture(opts)
    #arch = declarations_output_inhibit.Declare_Architecture(opts)
    #Warning! Make sure you know what declarations you are pulling from
    arch = declarations.Declare_Architecture(opts)
    #Build the forwards euler loop
    file_body_forwards = Forwards_Method_cupy_recovery_batch.Euler_Compiler(arch[0],arch[1],arch[2],opts,num_cells)
    #Compile a recovery-aware solve module in memory so this run path does not overwrite generated_solve_file.py.
    generated_solve_file = build_recovery_solver(file_body_forwards)


    ############
    #- Move the data loading to a seperate file and make it toggleable

    #Calculate the approximate spontaneous firing rate before fitting.
    # --- paths (MATLAB: cd(userpath); cd('../GitHub/.../OliverDataPlotting')) ---
    userpath = Path(r"C:\Users\ipboy\Documents\MATLAB")  # <-- change this
    plot_dir = userpath / "../GitHub/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting"
    plot_dir = plot_dir.resolve()
    os.chdir(plot_dir)

    # --- load MAT files ---
    mat = loadmat("all_units_info_with_polished_criteria_modified_perf.mat",variable_names=["all_data"],squeeze_me=True,struct_as_record=False)
    all_data = mat["all_data"]  # numpy array of MATLAB structs

    #cell_storage = []

    

    no_pre_holder = np.zeros((num_cells,10,29801))
    FRs = []

    #Get peak indicies
    peak_indices = []
    for n in range(num_cells):
        vals = []
        for m in range(4):
            vals2 = 0
            for k in range(10):
                if len(np.shape(all_data[n].ctrl_tar1_timestamps[k, m])) != 0:  #This just skips areas where there are no spikes
                    vals2 += np.shape(all_data[n].ctrl_tar1_timestamps[k, m])[0]
            vals.append(vals2)

        idx = np.where(vals == np.max(np.array(vals)))[0][0]
        vals = []
        peak_indices.append(idx)
       



    for k in range(num_cells):
        spike_times = all_data[k].ctrl_tar1_timestamps
        spike_times = spike_times[:, peak_indices[k]]  #Grab -90 degrees  #4.1.2026 -> Instead of grabbing -90 grab the indicy corresponding to the peak of the tuning curve.
        #cell_storage.append(spike_times)

        pre_zeros_list = []
        post_zeros_list = []

        for m in range(10):  # MATLAB: 1:10
            times = np.asarray(spike_times[m]).squeeze()
            no_pre_holder[k,m,np.round(times[(times < 2.98) & (times >= 0)]*(1000/opts['dt'])).astype(int)] = 1  # Convert to indices
            #pre_to_aprx_3s_holder[k,np.round(times[(times < 2.98)]*(1000/opts['dt'])+(1000/opts['dt'])).astype(int)] = 1  # Convert to indices
            #no_pre_holder[k,np.round(times[times < 0]*(1000/opts['dt'])).astype(int)] = 1  # Convert to indices

            pre_zeros_list.append(times[times < 0])
            #post_zeros_list.append(times[times >= 0])

        pre_zeros = np.concatenate(pre_zeros_list) if pre_zeros_list else np.array([])
        FR = pre_zeros.size / 10  # MATLAB: length(pre_zeros)/10
        FRs.append(FR)

    #n = 7  # MATLAB is 1-based
    #unit = all_data[n - 1]

    #spike_times = unit.ctrl_tar1_timestamps

    # Ensure we have "(:,1)" behavior: take first column if 2D
    #if isinstance(spike_times, np.ndarray) and spike_times.ndim == 2:
    #    spike_times = spike_times[:, 0]

    #pre_to_aprx_3s_holder = np.zeros((10,39801))


    
    #pre_to_aprx_3s_holder = np.zeros((10,29801))
    

    # pre_zeros_list = []
    # post_zeros_list = []
    # for k in range(10):  # MATLAB: 1:10
    #     times = np.asarray(spike_times[k]).squeeze()
    #     pre_to_aprx_3s_holder[k,np.round(times[(times < 2.98) & (times >= 0)]*(1000/opts['dt'])).astype(int)] = 1  # Convert to indices
    #     #pre_to_aprx_3s_holder[k,np.round(times[(times < 2.98)]*(1000/opts['dt'])+(1000/opts['dt'])).astype(int)] = 1  # Convert to indices
    #     #no_pre_holder[k,np.round(times[times < 0]*(1000/opts['dt'])).astype(int)] = 1  # Convert to indices

    #     pre_zeros_list.append(times[times < 0])
    #     post_zeros_list.append(times[times >= 0])
 
    # pre_zeros = np.concatenate(pre_zeros_list) if pre_zeros_list else np.array([])
    # FR = pre_zeros.size / 10  # MATLAB: length(pre_zeros)/10

    # print(FR)

    # -- Load in data
    #filename = f"C:/Users/ipboy/Documents/Github/ModelingEffort/Multi-Channel/Plotting/OliverDataPlotting/PicturesToFit/picture_fit{n}contra.mat"
    #data = loadmat(filename)['picture'].astype(np.float32)  #trials,timecourse
    
    #data = data[:,:,None]
    data = no_pre_holder[:,:,:,None]
    isi_bin_edges = make_isi_score_bin_edges(opts["sim_len"])
    data_isi_reference = precompute_data_isi_reference(no_pre_holder, isi_bin_edges)

    num_params = 8
    batch_size = opts['N_batch']
    p, fixed_params_per_cell = load_fixed_params_from_mat(
        START_PARAMS_PATH,
        batch_size,
        num_params,
        num_cells,
    )
    print(f"Loaded fixed non-refractory params from {START_PARAMS_PATH}")


    #spks = call_inputs(num_cells,FRs,batch_size)
    #on_spks = np.transpose(spks[f'locs_masker_None_target_0_on'][f'stimulus_0_poisson_spks'],(2,0,1))
    #off_spks = np.transpose(spks[f'locs_masker_None_target_0_off'][f'stimulus_0_poisson_spks'],(2,0,1))
    #noise = np.transpose(spks['noise_masker_None_target_0'],(0,1,4,2,3))

    #Repeat along the first axis to match the data shape #Note! You might be able to save some space by just using this across all tials... all trials are giong to have the bigger size anyways through so it might not matter much.
    #On second thought this should broadcast correctly.
    #on_spks= np.repeat(on_spks[None, ...], 220, axis=0) 
    #off_spks= np.repeat(off_spks[None, ...], 220, axis=0) 
    #noise= np.repeat(noise[None, ...], 220, axis=0) 

    #print('p1')
    #print(p[13:17,1:5])

    # Saving this full `noise` tensor to a v5 MAT-file will often exceed MATLAB's ~2GB limit.
    # Use the Python raster plot below instead of exporting the full array.

    #if plot_noise_toggle == 1:

    plot_noise_toggle = 1
    plot_noise_cell = 9      # MATLAB-style indexing (cell 1..num_cells)
    plot_noise_batch = 2     # MATLAB-style indexing (batch 1..N_batch)
    plot_noise_channel = 1   # MATLAB-style indexing (channel 1..chans)


    #from plot_noise_raster import plot_noise_raster
    #plot_noise_raster(noise,cell=plot_noise_cell,batch=plot_noise_batch,channel=plot_noise_channel,dt_ms=opts['dt'],one_indexed=True,)
    # Keep Adam-shaped arrays in the checkpoint for compatibility with the other runners.
    m = np.zeros((num_params,num_cells,batch_size))
    v = np.zeros((num_params,num_cells,batch_size))
    t = 0
    

    losses = []
    param_tracker = []
    start_epoch = 0

    best_output_per_cell = np.zeros((num_cells, opts["N_trials"], opts["N_channels"], opts["sim_len"]), dtype=np.int8)
    best_loss_per_cell = np.full((num_cells,), np.inf, dtype=np.float32)
    best_batch_id_per_cell = np.zeros((num_cells,), dtype=np.int32)
    best_epoch_per_cell = np.zeros((num_cells,), dtype=np.int32)
    best_params_per_cell = fixed_params_per_cell.copy()
    rng_refractory = np.random.default_rng(REFRACTORY_RANDOM_SEED)
    start_refractory_params = load_start_refractory_params_from_mat(START_PARAMS_PATH, num_cells)
    if start_refractory_params is None:
        refractory_param_center = estimate_initial_refractory_params(no_pre_holder)
        print("No refractory params found in the starting file; initialized refractory search from data ISIs.")
    else:
        refractory_param_center = start_refractory_params
        print("Initialized refractory search from best_refractory_params_per_cell in the starting file.")
    refractory_param_sigma = REFRACTORY_INITIAL_SIGMA.copy()
    initial_abs_steps, initial_tau_steps, initial_shape, initial_ceiling = refractory_params_to_values(refractory_param_center)
    best_refractory_params_per_cell = refractory_param_center.copy()
    best_refractory_abs_steps_per_cell = initial_abs_steps.astype(np.float32, copy=True)
    best_refractory_tau_steps_per_cell = initial_tau_steps.astype(np.float32, copy=True)
    best_refractory_shape_per_cell = initial_shape.astype(np.float32, copy=True)
    best_refractory_ceiling_per_cell = initial_ceiling.astype(np.float32, copy=True)
    best_recovery_funcs_per_cell = build_parametric_recovery_functions(
        refractory_param_center[:, None, :],
        recovery_funcs.shape[1],
    )[:, 0, :]
    best_isi_hist_l1_per_cell = np.full((num_cells,), np.nan, dtype=np.float32)
    best_isi_median_penalty_per_cell = np.full((num_cells,), np.nan, dtype=np.float32)
    best_isi_count_penalty_per_cell = np.full((num_cells,), np.nan, dtype=np.float32)
    best_isi_cv_penalty_per_cell = np.full((num_cells,), np.nan, dtype=np.float32)
    refractory_param_history = np.zeros((TOTAL_EPOCHS, num_cells, batch_size, 4), dtype=np.float32)
    refractory_score_history = np.zeros((TOTAL_EPOCHS, num_cells, batch_size), dtype=np.float32)

    if RESUME_CHECKPOINT_PATH:
        expected_shape = p.shape
        p, m, v, t, loaded_epoch, checkpoint_best_state = load_training_checkpoint(
            checkpoint_path=RESUME_CHECKPOINT_PATH,
            expected_shape=expected_shape,
        )
        start_epoch = loaded_epoch
        if checkpoint_best_state:
            loaded_best_output_per_cell = np.asarray(
                checkpoint_best_state.get("best_output_per_cell", best_output_per_cell),
                dtype=np.int8,
            )
            if loaded_best_output_per_cell.ndim == 3:
                loaded_best_output_per_cell = loaded_best_output_per_cell[:, :, None, :]
            best_output_per_cell = loaded_best_output_per_cell
            best_loss_per_cell = np.asarray(
                checkpoint_best_state.get("best_loss_per_cell", best_loss_per_cell),
                dtype=np.float32,
            ).reshape(num_cells)
            best_batch_id_per_cell = np.asarray(
                checkpoint_best_state.get("best_batch_id_per_cell", best_batch_id_per_cell),
                dtype=np.int32,
            ).reshape(num_cells)
            best_epoch_per_cell = np.asarray(
                checkpoint_best_state.get("best_epoch_per_cell", best_epoch_per_cell),
                dtype=np.int32,
            ).reshape(num_cells)
            best_params_per_cell = np.asarray(
                checkpoint_best_state.get("best_params_per_cell", best_params_per_cell),
                dtype=np.float32,
            )
            fixed_params_per_cell = np.asarray(
                checkpoint_best_state.get("fixed_params_per_cell", fixed_params_per_cell),
                dtype=np.float32,
            )
            best_refractory_params_per_cell = np.asarray(
                checkpoint_best_state.get("best_refractory_params_per_cell", best_refractory_params_per_cell),
                dtype=np.float32,
            ).reshape(num_cells, 4)
            best_refractory_abs_steps_per_cell = np.asarray(
                checkpoint_best_state.get("best_refractory_abs_steps_per_cell", best_refractory_abs_steps_per_cell),
                dtype=np.float32,
            ).reshape(num_cells)
            best_refractory_tau_steps_per_cell = np.asarray(
                checkpoint_best_state.get("best_refractory_tau_steps_per_cell", best_refractory_tau_steps_per_cell),
                dtype=np.float32,
            ).reshape(num_cells)
            best_refractory_shape_per_cell = np.asarray(
                checkpoint_best_state.get("best_refractory_shape_per_cell", best_refractory_shape_per_cell),
                dtype=np.float32,
            ).reshape(num_cells)
            best_refractory_ceiling_per_cell = np.asarray(
                checkpoint_best_state.get("best_refractory_ceiling_per_cell", best_refractory_ceiling_per_cell),
                dtype=np.float32,
            ).reshape(num_cells)
            best_recovery_funcs_per_cell = np.asarray(
                checkpoint_best_state.get("best_recovery_funcs_per_cell", best_recovery_funcs_per_cell),
                dtype=np.float32,
            )
            best_isi_hist_l1_per_cell = np.asarray(
                checkpoint_best_state.get("best_isi_hist_l1_per_cell", best_isi_hist_l1_per_cell),
                dtype=np.float32,
            ).reshape(num_cells)
            best_isi_median_penalty_per_cell = np.asarray(
                checkpoint_best_state.get("best_isi_median_penalty_per_cell", best_isi_median_penalty_per_cell),
                dtype=np.float32,
            ).reshape(num_cells)
            best_isi_count_penalty_per_cell = np.asarray(
                checkpoint_best_state.get("best_isi_count_penalty_per_cell", best_isi_count_penalty_per_cell),
                dtype=np.float32,
            ).reshape(num_cells)
            best_isi_cv_penalty_per_cell = np.asarray(
                checkpoint_best_state.get("best_isi_cv_penalty_per_cell", best_isi_cv_penalty_per_cell),
                dtype=np.float32,
            ).reshape(num_cells)
            refractory_param_center = best_refractory_params_per_cell.copy()
        print(
            f"Resumed from checkpoint '{RESUME_CHECKPOINT_PATH}' "
            f"(saved after epoch {loaded_epoch})."
        )

    target_dict = call_strfs(p, batch_size, num_cells)
    spks = call_inputs(num_cells, FRs, batch_size, target_dict)

    on_spks = np.transpose(spks[f'locs_masker_None_target_0_on'][f'stimulus_0_poisson_spks'], (3, 4, 2, 1, 0))
    off_spks = np.transpose(spks[f'locs_masker_None_target_0_off'][f'stimulus_0_poisson_spks'], (3, 4, 2, 1, 0))
    rate_on = spks[f'locs_masker_None_target_0_on'][f'stimulus_0_rate']
    rate_off = spks[f'locs_masker_None_target_0_off'][f'stimulus_0_rate']
    rate_on_deriv = spks[f'locs_masker_None_target_0_on'][f'stimulus_0_rate_deriv']
    rate_off_deriv = spks[f'locs_masker_None_target_0_off'][f'stimulus_0_rate_deriv']
    noise = np.transpose(spks['noise_masker_None_target_0'], (0, 1, 4, 3, 2))

    if FIX_INPUTS_ACROSS_REFRACTORY_CANDIDATES:
        on_spks = repeat_first_batch(on_spks, batch_axis=1, batch_size=batch_size)
        off_spks = repeat_first_batch(off_spks, batch_axis=1, batch_size=batch_size)
        noise = repeat_first_batch(noise, batch_axis=1, batch_size=batch_size)
        rate_on = repeat_first_batch(rate_on, batch_axis=2, batch_size=batch_size)
        rate_off = repeat_first_batch(rate_off, batch_axis=2, batch_size=batch_size)
        rate_on_deriv = repeat_first_batch(rate_on_deriv, batch_axis=2, batch_size=batch_size)
        rate_off_deriv = repeat_first_batch(rate_off_deriv, batch_axis=2, batch_size=batch_size)
        print("Using identical input/noise streams across refractory candidates.")

    best_loss = 1e32

    best_output = []

    best_params = []

    start = time.perf_counter()

    for epoch in range(start_epoch, TOTAL_EPOCHS):

        print(range(start_epoch, TOTAL_EPOCHS))
        refractory_param_candidates = propose_refractory_params(
            refractory_param_center,
            refractory_param_sigma,
            batch_size,
            rng_refractory,
        )
        recovery_funcs_this_epoch = build_parametric_recovery_functions(
            refractory_param_candidates,
            recovery_funcs.shape[1],
        )
        refractory_param_history[epoch] = refractory_param_candidates

        output, grads, on_track, off_track, loss_holder = generated_solve_file.solve_run(on_spks,off_spks,noise,rate_on,rate_off,rate_on_deriv,rate_off_deriv,data,p,recovery_funcs_this_epoch) #Python Verison to build
        refractory_score, refractory_score_parts = score_refractory_candidates(
            output,
            data_isi_reference,
            isi_bin_edges,
        )
        losses.append(float(np.nanmedian(refractory_score)))
        param_tracker.append(fixed_params_per_cell.copy())

        print(
            f'ISI-only refractory objective median: {np.nanmedian(refractory_score):.3f} '
            f'---- Epoch: {epoch}'
        )
        refractory_score_history[epoch] = refractory_score
        best_batch_idx_this_epoch = np.argmin(refractory_score, axis=1)
        best_loss_this_epoch = refractory_score[np.arange(num_cells), best_batch_idx_this_epoch]
        improved_cells = best_loss_this_epoch < best_loss_per_cell
        best_params_this_epoch = refractory_param_candidates[
            np.arange(num_cells),
            best_batch_idx_this_epoch,
            :,
        ].astype(np.float32, copy=False)
        refractory_param_center = (
            REFRACTORY_CENTER_MOMENTUM * refractory_param_center
            + (1.0 - REFRACTORY_CENTER_MOMENTUM) * best_params_this_epoch
        ).astype(np.float32, copy=False)
        refractory_param_center = clip_refractory_params(refractory_param_center)
        refractory_param_sigma = np.maximum(
            refractory_param_sigma * REFRACTORY_SIGMA_DECAY,
            REFRACTORY_MIN_SIGMA,
        )

        if np.any(improved_cells):
            improved_idx = np.where(improved_cells)[0]
            chosen_batch_idx = best_batch_idx_this_epoch[improved_idx]
            chosen_refractory_params = refractory_param_candidates[improved_idx, chosen_batch_idx, :]
            chosen_abs, chosen_tau, chosen_shape, chosen_ceiling = refractory_params_to_values(chosen_refractory_params)

            best_output_per_cell[improved_idx, :, :, :] = np.asarray(
                output[improved_idx, chosen_batch_idx, :, :, :],
                dtype=np.int8,
            )
            best_loss_per_cell[improved_idx] = best_loss_this_epoch[improved_idx]
            best_batch_id_per_cell[improved_idx] = chosen_batch_idx.astype(np.int32) + 1
            best_epoch_per_cell[improved_idx] = epoch + 1
            best_params_per_cell[:, improved_idx] = fixed_params_per_cell[:, improved_idx].astype(np.float32)
            best_refractory_params_per_cell[improved_idx] = chosen_refractory_params
            best_refractory_abs_steps_per_cell[improved_idx] = chosen_abs
            best_refractory_tau_steps_per_cell[improved_idx] = chosen_tau
            best_refractory_shape_per_cell[improved_idx] = chosen_shape
            best_refractory_ceiling_per_cell[improved_idx] = chosen_ceiling
            best_recovery_funcs_per_cell[improved_idx] = recovery_funcs_this_epoch[improved_idx, chosen_batch_idx]
            best_isi_hist_l1_per_cell[improved_idx] = refractory_score_parts["isi_hist_l1"][improved_idx, chosen_batch_idx]
            best_isi_median_penalty_per_cell[improved_idx] = refractory_score_parts["isi_median_penalty"][improved_idx, chosen_batch_idx]
            best_isi_count_penalty_per_cell[improved_idx] = refractory_score_parts["isi_count_penalty"][improved_idx, chosen_batch_idx]
            best_isi_cv_penalty_per_cell[improved_idx] = refractory_score_parts["isi_cv_penalty"][improved_idx, chosen_batch_idx]

            print(
                f'Updated best output for {len(improved_idx)} cells '
                f'(median best refractory objective: {np.nanmedian(best_loss_per_cell):.3f}, '
                f'median abs/tau/shape/ceil: '
                f'{np.nanmedian(best_refractory_abs_steps_per_cell):.1f}/'
                f'{np.nanmedian(best_refractory_tau_steps_per_cell):.1f}/'
                f'{np.nanmedian(best_refractory_shape_per_cell):.2f}/'
                f'{np.nanmedian(best_refractory_ceiling_per_cell):.2f})'
            )

        epoch_one_indexed = epoch + 1
        if epoch_one_indexed % CHECKPOINT_EVERY_EPOCHS == 0:
            save_training_checkpoint(
                epoch_one_indexed=epoch_one_indexed,
                p=p,
                m=m,
                v=v,
                t=t,
                losses=losses,
                best_output_per_cell=best_output_per_cell,
                best_loss_per_cell=best_loss_per_cell,
                best_batch_id_per_cell=best_batch_id_per_cell,
                best_epoch_per_cell=best_epoch_per_cell,
                best_params_per_cell=best_params_per_cell,
                best_refractory_params_per_cell=best_refractory_params_per_cell,
                best_refractory_abs_steps_per_cell=best_refractory_abs_steps_per_cell,
                best_refractory_tau_steps_per_cell=best_refractory_tau_steps_per_cell,
                best_refractory_shape_per_cell=best_refractory_shape_per_cell,
                best_refractory_ceiling_per_cell=best_refractory_ceiling_per_cell,
                best_recovery_funcs_per_cell=best_recovery_funcs_per_cell,
                refractory_param_history=refractory_param_history[:epoch_one_indexed],
                refractory_score_history=refractory_score_history[:epoch_one_indexed],
                best_isi_hist_l1_per_cell=best_isi_hist_l1_per_cell,
                best_isi_median_penalty_per_cell=best_isi_median_penalty_per_cell,
                best_isi_count_penalty_per_cell=best_isi_count_penalty_per_cell,
                best_isi_cv_penalty_per_cell=best_isi_cv_penalty_per_cell,
                fixed_params_per_cell=fixed_params_per_cell,
            )

        #print('p3')
        #print(p[13:17,1:5])

        #print(p[0][0])

        #print(np.shape(loss))
        #print(np.shape(p))

        


    elapsed = time.perf_counter() - start
    print(f"{elapsed:.2f} s")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    #savemat(f"output_compressed_Eprop_{timestamp}.mat", {"output": output, "losses":losses, "params" : param_tracker,  "best_loss" : np.asarray(best_loss, dtype=np.float32),"best_output" : np.asarray(best_output, dtype=np.float32),"best_params" : np.asarray(best_params, dtype=np.float32)}, do_compression=True)
    savemat(
        f"output_compressed_Eprop_All_cells_recovery_fit_refractory_only_{timestamp}.mat",
        {
            "output": output,
            "losses": losses,
            "params": p,
            "fixed_params_per_cell": fixed_params_per_cell,
            "best_output_per_cell": best_output_per_cell,
            "best_loss_per_cell": best_loss_per_cell,
            "best_batch_id_per_cell": best_batch_id_per_cell,
            "best_epoch_per_cell": best_epoch_per_cell,
            "best_params_per_cell": best_params_per_cell,
            "best_refractory_params_per_cell": best_refractory_params_per_cell,
            "best_refractory_abs_steps_per_cell": best_refractory_abs_steps_per_cell,
            "best_refractory_tau_steps_per_cell": best_refractory_tau_steps_per_cell,
            "best_refractory_shape_per_cell": best_refractory_shape_per_cell,
            "best_refractory_ceiling_per_cell": best_refractory_ceiling_per_cell,
            "best_recovery_funcs_per_cell": best_recovery_funcs_per_cell,
            "best_isi_hist_l1_per_cell": best_isi_hist_l1_per_cell,
            "best_isi_median_penalty_per_cell": best_isi_median_penalty_per_cell,
            "best_isi_count_penalty_per_cell": best_isi_count_penalty_per_cell,
            "best_isi_cv_penalty_per_cell": best_isi_cv_penalty_per_cell,
            "refractory_param_history": refractory_param_history,
            "refractory_score_history": refractory_score_history,
            "recovery_template_from_file": recovery_funcs,
            "isi_score_bin_edges": isi_bin_edges,
            "start_params_path": str(START_PARAMS_PATH),
            "isi_objective_note": "Refractory-only fit. Fixed model params; no PSTH loss or Adam update used.",
        },
        do_compression=True,
    )

    # # ============== PARAMETER EVOLUTION PLOT ==============
    # # param_tracker is a list of arrays, each with shape (1, 100) or (num_params, batch_size)
    # # Stack into array: (num_epochs, num_params, batch_size)
    # param_array = np.array(param_tracker)
    # print(f"param_array shape: {param_array.shape}")
    
    # fig_params, axes = plt.subplots(2, 2, figsize=(14, 10))
    # fig_params.suptitle('Parameter Evolution Across Epochs (All 100 Batch Members)', fontsize=14)
    
    # # Plot 1: Each parameter trajectory over epochs (all 100 batch members)
    # ax = axes[0, 0]
    # epochs_x = np.arange(len(param_tracker))
    # for batch_idx in range(param_array.shape[2]):
    #     ax.plot(epochs_x, param_array[:, 0, batch_idx], alpha=0.3, linewidth=1)
    # ax.set_xlabel('Epoch')
    # ax.set_ylabel('Parameter Value')
    # ax.set_title('All 100 Parameter Trajectories')
    # ax.grid(True, alpha=0.3)
    
    # # Plot 2: Parameter distribution at each epoch (box plot style)
    # # ax = axes[0, 1]
    # # param_at_epochs = [param_array[e, 0, :] for e in range(len(param_tracker))]
    # # bp = ax.boxplot(param_at_epochs, positions=epochs_x)
    # # ax.set_xlabel('Epoch')
    # # ax.set_ylabel('Parameter Value')
    # # ax.set_title('Parameter Distribution per Epoch')
    # # ax.grid(True, alpha=0.3)
    
    # # Plot 3: Parameter changes (delta) per epoch
    # ax = axes[1, 0]
    # if len(param_tracker) > 1:
    #     for batch_idx in range(param_array.shape[2]):
    #         deltas = np.diff(param_array[:, 0, batch_idx])
    #         ax.plot(epochs_x[1:], deltas, alpha=0.3, linewidth=1)
    #     ax.axhline(y=0, color='red', linestyle='--', linewidth=1, label='Zero change')
    #     ax.set_xlabel('Epoch')
    #     ax.set_ylabel('Parameter Change (Δp)')
    #     ax.set_title('Parameter Step Size per Epoch')
    #     ax.legend()
    # ax.grid(True, alpha=0.3)
    
    # # Plot 4: Initial vs Final parameter values (scatter)
    # ax = axes[1, 1]
    # initial_params = param_array[0, 0, :]
    # final_params = param_array[-1, 0, :]
    # ax.scatter(initial_params, final_params, alpha=0.5, s=30)
    # # Add diagonal line for reference
    # min_val = min(initial_params.min(), final_params.min())
    # max_val = max(initial_params.max(), final_params.max())
    # ax.plot([min_val, max_val], [min_val, max_val], 'r--', label='No change line')
    # ax.set_xlabel('Initial Parameter Value')
    # ax.set_ylabel('Final Parameter Value')
    # ax.set_title('Initial vs Final Parameters')
    # ax.legend()
    # ax.grid(True, alpha=0.3)
    
    # plt.tight_layout()
    # plt.savefig('parameter_evolution.png', dpi=150, bbox_inches='tight')
    # plt.show()
    # print("Saved parameter evolution plot to parameter_evolution.png")

if __name__ == "__main__":


    run_sim = runSimulation()
