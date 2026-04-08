from __future__ import annotations

import base64
from io import BytesIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import to_rgb
from scipy.io import loadmat
from scipy.io.matlab import mat_struct
from sklearn.preprocessing import StandardScaler

try:
    from umap import UMAP
except ImportError as exc:
    raise ImportError(
        "UMAP could not be imported. Install `umap-learn` in the Python environment "
        "you use to run this script."
    ) from exc


sns.set_theme(style="whitegrid", context="notebook")

OUTPUT_MAT_PATH = Path(
    r"C:\Users\ipboy\Documents\GitHub\ModelingEffort\Multi-Channel\Plotting\OliverDataPlotting\output_compressed_Eprop_All_cells_20260402_090929.mat"
)
ALL_UNITS_PATH = Path(
    r"C:\Users\ipboy\Documents\GitHub\ModelingEffort\Multi-Channel\Plotting\OliverDataPlotting\all_units_info_with_polished_criteria_modified_perf.mat"
)

OUTPUT_DIR = Path(__file__).resolve().parent / "umap_outputs"

DT_MS = 0.1
TIMESTEPS = 29801
PREVIEW_MAX_TIME_BINS = 400
PREVIEW_ROW_SCALE = 10
UMAP_RANDOM_STATE = 42
UMAP_N_NEIGHBORS = 5
UMAP_MIN_DIST = 1e-6
MAKE_BOKEH_HOVER = True
TOPOGRAPHY_GRID_SIZE = 220
TOPOGRAPHY_ALPHA = 0.38

LAYER_ORDER = ["L2/3", "L4", "L5/6", "NaN"]
LAYER_PALETTE = {
    "L2/3": "#396ab1",
    "L4": "#3e9651",
    "L5/6": "#cc2529",
    "NaN": "#a0a0a0",
}
FALLBACK_LAYER_COLOR = "#7f7f7f"

# Fill these in if you want manuscript-style names instead of param_1, param_2, ...
PARAMETER_NAME_OVERRIDES: dict[int, str] = {}


def load_all_data(all_units_path: Path):
    mat = loadmat(
        all_units_path,
        variable_names=["all_data"],
        squeeze_me=True,
        struct_as_record=False,
    )
    if "all_data" not in mat:
        keys = [k for k in mat.keys() if not k.startswith("__")]
        raise KeyError(
            f"Expected 'all_data' in {all_units_path}. Available keys: {keys}"
        )
    return mat["all_data"]


def load_output_bundle(output_mat_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mat = loadmat(output_mat_path, variable_names=["output", "losses", "params"])
    if "output" not in mat or "losses" not in mat or "params" not in mat:
        keys = [k for k in mat.keys() if not k.startswith("__")]
        raise KeyError(
            f"Expected 'output', 'losses', and 'params' in {output_mat_path}. Available keys: {keys}"
        )
    output = np.asarray(mat["output"])
    losses = np.asarray(mat["losses"])
    params = np.asarray(mat["params"])
    return output, losses, params


def best_batch_indices_from_final_epoch(losses: np.ndarray) -> np.ndarray:
    losses = np.asarray(losses)
    if losses.ndim != 4:
        raise ValueError(
            "Expected losses to have shape (epoch, metric, cell, batch), "
            f"got {losses.shape}"
        )

    final_epoch_binned_sse = losses[-1, 1, :, :]
    return np.argmin(final_epoch_binned_sse, axis=-1) + 1


def build_parameter_names(param_count: int) -> list[str]:
    return [
        PARAMETER_NAME_OVERRIDES.get(index, f"param_{index + 1}")
        for index in range(param_count)
    ]


def normalize_layer_label(value: object) -> str:
    if value is None:
        return "NaN"

    label = str(value).strip()
    if label == "" or label.lower() == "nan":
        return "NaN"
    return label


def ordered_layers(layer_values: pd.Series) -> list[str]:
    unique_layers = list(pd.unique(layer_values.astype(str)))
    known_layers = [layer for layer in LAYER_ORDER if layer in unique_layers]
    extra_layers = sorted(layer for layer in unique_layers if layer not in known_layers)
    return known_layers + extra_layers


def palette_for_layers(layer_labels: list[str]) -> dict[str, str]:
    return {
        layer: LAYER_PALETTE.get(layer, FALLBACK_LAYER_COLOR)
        for layer in layer_labels
    }


def collect_best_fit_parameters(
    params: np.ndarray,
    best_batch_indices_1based: np.ndarray,
) -> np.ndarray:
    params = np.asarray(params)
    if params.ndim != 4:
        raise ValueError(
            "Expected params to have shape (epoch, parameter, cell, batch), "
            f"got {params.shape}"
        )

    epoch_idx = params.shape[0] - 1
    unit_count = min(params.shape[2], best_batch_indices_1based.size)

    best_params = []
    for unit_idx in range(unit_count):
        batch_idx = int(best_batch_indices_1based[unit_idx]) - 1
        if batch_idx < 0 or batch_idx >= params.shape[3]:
            raise IndexError(
                f"Best batch index {batch_idx + 1} is out of bounds for unit {unit_idx + 1}"
            )
        best_params.append(np.asarray(params[epoch_idx, :, unit_idx, batch_idx], dtype=float))

    return np.vstack(best_params)


def extract_layer_labels(all_data, unit_count: int) -> list[str]:
    labels: list[str] = []
    for unit_idx in range(unit_count):
        unit = all_data[unit_idx]
        if not isinstance(unit, mat_struct):
            raise TypeError(
                f"Unexpected MATLAB struct format for unit {unit_idx + 1}: {type(unit)}"
            )
        labels.append(normalize_layer_label(getattr(unit, "layer", "NaN")))
    return labels


def build_analysis_dataframe(
    all_data,
    best_batch_indices_1based: np.ndarray,
    params: np.ndarray,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    best_params = collect_best_fit_parameters(params, best_batch_indices_1based)
    unit_count = min(len(all_data), best_params.shape[0], best_batch_indices_1based.size)

    if unit_count < best_params.shape[0]:
        best_params = best_params[:unit_count, :]
        best_batch_indices_1based = best_batch_indices_1based[:unit_count]

    param_cols = build_parameter_names(best_params.shape[1])
    df = pd.DataFrame(best_params, columns=param_cols)
    df.insert(0, "unit_id", np.arange(1, unit_count + 1, dtype=int))
    df["best_batch"] = best_batch_indices_1based.astype(int)
    df["layer"] = extract_layer_labels(all_data, unit_count)

    finite_mask = np.isfinite(df[param_cols].to_numpy()).all(axis=1)
    dropped = int((~finite_mask).sum())
    if dropped:
        print(f"Dropping {dropped} unit(s) with non-finite parameter values before plotting.")
        df = df.loc[finite_mask].reset_index(drop=True)

    layer_labels = ordered_layers(df["layer"])
    return df, param_cols, layer_labels


def compute_umap_embedding(
    df: pd.DataFrame,
    param_cols: list[str],
) -> pd.DataFrame:
    if len(df) < 2:
        raise ValueError("Need at least two cells to compute a UMAP embedding.")

    scaled = StandardScaler().fit_transform(df[param_cols].to_numpy())
    n_neighbors = max(2, min(UMAP_N_NEIGHBORS, len(df) - 1))
    reducer = UMAP(
        n_neighbors=n_neighbors,
        min_dist=UMAP_MIN_DIST,
        random_state=UMAP_RANDOM_STATE,
    )
    embedding = reducer.fit_transform(scaled)

    embedded_df = df.copy()
    embedded_df["umap_1"] = embedding[:, 0]
    embedded_df["umap_2"] = embedding[:, 1]
    return embedded_df


def ensure_trial_time(arr: np.ndarray) -> np.ndarray:
    arr = np.asarray(arr)
    arr = np.squeeze(arr)

    if arr.ndim == 1:
        return arr[None, :]

    if arr.ndim == 2:
        if arr.shape[0] == 10:
            return arr
        if arr.shape[1] == 10:
            return arr.T
        return arr

    candidate_axes = [axis for axis, size in enumerate(arr.shape) if size == 10]
    if not candidate_axes:
        raise ValueError(f"Could not find a trial axis of length 10 in shape {arr.shape}")

    moved = np.moveaxis(arr, candidate_axes[0], 0)
    moved = np.squeeze(moved)
    if moved.ndim == 1:
        moved = moved[None, :]
    if moved.ndim != 2:
        raise ValueError(f"Could not coerce shape {arr.shape} to (trials, time)")
    return moved


def extract_numeric_spike_times(value: object) -> np.ndarray:
    arr = np.asarray(value)
    arr = np.squeeze(arr)

    if arr.dtype != object:
        return np.atleast_1d(arr).astype(float, copy=False).ravel()

    if arr.ndim == 0:
        return extract_numeric_spike_times(arr.item())

    if arr.size == 0:
        return np.empty(0, dtype=float)

    flattened_parts = [
        extract_numeric_spike_times(item)
        for item in arr.reshape(-1)
    ]
    non_empty_parts = [part for part in flattened_parts if part.size > 0]
    if not non_empty_parts:
        return np.empty(0, dtype=float)

    return np.concatenate(non_empty_parts).astype(float, copy=False)


def compute_peak_direction_indices(all_data) -> list[int]:
    peak_indices: list[int] = []

    for unit in all_data:
        if not isinstance(unit, mat_struct):
            raise TypeError("Unexpected MATLAB struct format in all_data")

        vals = []
        for direction_idx in range(4):
            count = 0
            for trial_idx in range(10):
                times = extract_numeric_spike_times(
                    unit.ctrl_tar1_timestamps[trial_idx, direction_idx]
                )
                count += times.size
            vals.append(count)

        peak_indices.append(int(np.argmax(vals)))

    return peak_indices


def build_target_raster(
    all_data,
    unit_index_1based: int,
    peak_indices: list[int],
    timesteps: int,
    dt_ms: float,
) -> np.ndarray:
    unit = all_data[unit_index_1based - 1]
    if not isinstance(unit, mat_struct):
        raise TypeError("Unexpected MATLAB struct format for selected unit")

    direction_idx = peak_indices[unit_index_1based - 1]

    raster = np.zeros((10, timesteps), dtype=np.float64)
    duration_s = timesteps * dt_ms / 1000.0

    for trial_idx in range(10):
        times = extract_numeric_spike_times(
            unit.ctrl_tar1_timestamps[trial_idx, direction_idx]
        )
        times = times[(times >= 0) & (times < duration_s)]
        inds = np.rint(times * (1000.0 / dt_ms)).astype(int)
        inds = inds[(inds >= 0) & (inds < timesteps)]
        raster[trial_idx, inds] = 1.0

    return raster


def load_model_raster(
    output: np.ndarray,
    unit_index_1based: int,
    batch_index_1based: int,
) -> np.ndarray:
    unit_idx = unit_index_1based - 1
    batch_idx = batch_index_1based - 1

    if output.ndim == 5:
        selected = output[unit_idx, batch_idx, :, :, :]
    elif output.ndim == 4:
        selected = output[unit_idx, :, :, :]
    else:
        raise ValueError(f"Unsupported output shape: {output.shape}")

    return ensure_trial_time(selected).astype(np.float64)


def downsample_raster_for_preview(
    raster: np.ndarray,
    max_time_bins: int = PREVIEW_MAX_TIME_BINS,
) -> np.ndarray:
    raster = np.asarray(raster, dtype=float)
    if raster.ndim != 2:
        raise ValueError(f"Expected raster shape (trials, time), got {raster.shape}")

    if raster.shape[1] <= max_time_bins:
        preview = raster
    else:
        step = int(np.ceil(raster.shape[1] / max_time_bins))
        pooled = []
        for start in range(0, raster.shape[1], step):
            stop = min(start + step, raster.shape[1])
            pooled.append(raster[:, start:stop].max(axis=1))
        preview = np.stack(pooled, axis=1)

    return np.repeat(preview, PREVIEW_ROW_SCALE, axis=0)


def raster_to_data_uri(raster: np.ndarray) -> str:
    preview = downsample_raster_for_preview(raster)
    buffer = BytesIO()
    plt.imsave(buffer, preview, cmap="binary", format="png", origin="lower", vmin=0.0, vmax=1.0)
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def build_hover_raster_previews(
    all_data,
    unit_ids_1based: np.ndarray,
    peak_indices: list[int],
) -> list[str]:
    previews: list[str] = []
    for unit_id in unit_ids_1based:
        try:
            raster = build_target_raster(
                all_data=all_data,
                unit_index_1based=int(unit_id),
                peak_indices=peak_indices,
                timesteps=TIMESTEPS,
                dt_ms=DT_MS,
            )
            previews.append(raster_to_data_uri(raster))
        except Exception as exc:
            print(f"Warning: could not build raster preview for unit {int(unit_id)}: {exc}")
            previews.append(raster_to_data_uri(np.zeros((10, 10), dtype=float)))
    return previews


def build_model_hover_raster_previews(
    output: np.ndarray,
    unit_ids_1based: np.ndarray,
    best_batch_indices_1based: np.ndarray,
) -> list[str]:
    previews: list[str] = []
    for unit_id, batch_index in zip(unit_ids_1based, best_batch_indices_1based):
        try:
            raster = load_model_raster(
                output=output,
                unit_index_1based=int(unit_id),
                batch_index_1based=int(batch_index),
            )
            previews.append(raster_to_data_uri(raster))
        except Exception as exc:
            print(
                "Warning: could not build model raster preview for "
                f"unit {int(unit_id)} batch {int(batch_index)}: {exc}"
            )
            previews.append(raster_to_data_uri(np.zeros((10, 10), dtype=float)))
    return previews


def compute_plot_bounds(
    df: pd.DataFrame,
    padding_fraction: float = 0.08,
) -> tuple[float, float, float, float]:
    x_min = float(df["umap_1"].min())
    x_max = float(df["umap_1"].max())
    y_min = float(df["umap_2"].min())
    y_max = float(df["umap_2"].max())

    x_span = max(x_max - x_min, 1e-6)
    y_span = max(y_max - y_min, 1e-6)
    x_pad = x_span * padding_fraction
    y_pad = y_span * padding_fraction

    return x_min - x_pad, x_max + x_pad, y_min - y_pad, y_max + y_pad


def estimate_topography_bandwidth(x_values: np.ndarray, y_values: np.ndarray) -> float:
    points = np.column_stack([x_values, y_values]).astype(float)
    if len(points) < 2:
        return 1.0

    diffs = points[:, None, :] - points[None, :, :]
    dists = np.sqrt(np.sum(diffs ** 2, axis=2))
    np.fill_diagonal(dists, np.inf)
    nearest_neighbor = np.min(dists, axis=1)
    finite_nn = nearest_neighbor[np.isfinite(nearest_neighbor)]

    span = max(float(np.ptp(x_values)), float(np.ptp(y_values)), 1e-6)
    if finite_nn.size == 0:
        return max(span / 12.0, 1e-3)

    return max(float(np.median(finite_nn)) * 2.0, span / 18.0, 1e-3)


def build_layer_topography_rgba(
    df: pd.DataFrame,
    layer_labels: list[str],
) -> tuple[np.ndarray, float, float, float, float]:
    palette = palette_for_layers(layer_labels)
    x_values = df["umap_1"].to_numpy(dtype=float)
    y_values = df["umap_2"].to_numpy(dtype=float)
    colors = np.vstack([np.asarray(to_rgb(palette[layer]), dtype=float) for layer in df["layer"]])

    x_min, x_max, y_min, y_max = compute_plot_bounds(df)
    xs = np.linspace(x_min, x_max, TOPOGRAPHY_GRID_SIZE)
    ys = np.linspace(y_min, y_max, TOPOGRAPHY_GRID_SIZE)
    grid_x, grid_y = np.meshgrid(xs, ys)

    bandwidth = estimate_topography_bandwidth(x_values, y_values)
    bandwidth_sq = bandwidth ** 2

    density = np.zeros_like(grid_x, dtype=float)
    rgb_numerator = np.zeros(grid_x.shape + (3,), dtype=float)

    for x0, y0, rgb in zip(x_values, y_values, colors):
        weight = np.exp(-0.5 * (((grid_x - x0) ** 2 + (grid_y - y0) ** 2) / bandwidth_sq))
        density += weight
        rgb_numerator += weight[..., None] * rgb

    safe_density = np.maximum(density, 1e-12)
    rgb = rgb_numerator / safe_density[..., None]

    positive_density = density[density > 0]
    if positive_density.size:
        reference = float(np.quantile(positive_density, 0.98))
    else:
        reference = 1.0
    alpha = np.clip(density / max(reference, 1e-12), 0.0, 1.0) ** 0.85
    alpha *= TOPOGRAPHY_ALPHA

    rgba = np.zeros(grid_x.shape + (4,), dtype=np.uint8)
    rgba[..., :3] = np.clip(rgb * 255.0, 0.0, 255.0).astype(np.uint8)
    rgba[..., 3] = np.clip(alpha * 255.0, 0.0, 255.0).astype(np.uint8)

    rgba_uint32 = np.ascontiguousarray(rgba).view(dtype=np.uint32).reshape(rgba.shape[:2])
    return rgba_uint32, x_min, y_min, x_max - x_min, y_max - y_min


def save_parameter_pairplot(
    df: pd.DataFrame,
    param_cols: list[str],
    layer_labels: list[str],
    save_path: Path,
) -> None:
    palette = palette_for_layers(layer_labels)

    if len(param_cols) == 1:
        fig, ax = plt.subplots(figsize=(8, 5))
        sns.histplot(
            data=df,
            x=param_cols[0],
            hue="layer",
            hue_order=layer_labels,
            palette=palette,
            multiple="stack",
            ax=ax,
        )
        ax.set_title("Best-fit parameter distribution by cortical layer")
        fig.tight_layout()
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        return

    grid = sns.pairplot(
        df,
        vars=param_cols,
        hue="layer",
        hue_order=layer_labels,
        palette=palette,
        corner=True,
        diag_kind="hist",
        plot_kws={"s": 48, "alpha": 0.85, "edgecolor": "white", "linewidth": 0.4},
        diag_kws={"alpha": 0.7},
    )
    grid.figure.suptitle("Best-fit parameter pairplot colored by cortical layer", y=1.02)
    grid.figure.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(grid.figure)


def save_seaborn_embedding_plot(
    df: pd.DataFrame,
    layer_labels: list[str],
    save_path: Path,
) -> None:
    palette = palette_for_layers(layer_labels)
    fig, ax = plt.subplots(figsize=(11, 9))
    sns.scatterplot(
        data=df,
        x="umap_1",
        y="umap_2",
        hue="layer",
        hue_order=layer_labels,
        palette=palette,
        s=70,
        edgecolor="white",
        linewidth=0.5,
        alpha=0.9,
        ax=ax,
    )
    ax.set_title("UMAP projection of best-fit parameters colored by cortical layer")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(title="Layer", frameon=True)
    fig.tight_layout()
    fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_bokeh_hover_embedding(
    df: pd.DataFrame,
    all_data,
    peak_indices: list[int],
    layer_labels: list[str],
    save_path: Path,
) -> None:
    try:
        from bokeh.io import output_file, save
        from bokeh.models import CategoricalColorMapper, ColumnDataSource, HoverTool
        from bokeh.plotting import figure
    except ImportError as exc:
        print(f"Skipping interactive hover plot because Bokeh is unavailable: {exc}")
        return

    palette = palette_for_layers(layer_labels)
    hover_df = df.copy()
    hover_df["raster_png"] = build_hover_raster_previews(
        all_data=all_data,
        unit_ids_1based=hover_df["unit_id"].to_numpy(dtype=int),
        peak_indices=peak_indices,
    )

    source = ColumnDataSource(hover_df)
    color_mapper = CategoricalColorMapper(
        factors=layer_labels,
        palette=[palette[layer] for layer in layer_labels],
    )

    hover_html = """
    <div style="width: 360px;">
      <div><strong>Unit:</strong> @unit_id</div>
      <div><strong>Layer:</strong> @layer</div>
      <div><strong>Best batch:</strong> @best_batch</div>
      <div style="margin-top: 8px;">
        <img
          src="@raster_png"
          alt="Raster preview"
          style="width: 320px; border: 1px solid #d0d0d0;"
        />
      </div>
    </div>
    """

    plot = figure(
        title="UMAP projection with target raster previews on hover",
        width=950,
        height=760,
        tools="pan,wheel_zoom,box_zoom,reset,save,hover",
        active_scroll="wheel_zoom",
    )
    plot.circle(
        x="umap_1",
        y="umap_2",
        source=source,
        size=9,
        alpha=0.9,
        line_color="white",
        line_width=0.5,
        fill_color={"field": "layer", "transform": color_mapper},
        legend_field="layer",
    )
    plot.xaxis.axis_label = "UMAP 1"
    plot.yaxis.axis_label = "UMAP 2"
    plot.legend.location = "top_right"
    plot.legend.click_policy = "hide"

    hover = plot.select_one(HoverTool)
    hover.tooltips = hover_html
    hover.point_policy = "follow_mouse"

    output_file(save_path, title="UMAP layer hover")
    save(plot)


def save_bokeh_hover_embedding_with_model_and_topography(
    df: pd.DataFrame,
    all_data,
    output: np.ndarray,
    peak_indices: list[int],
    layer_labels: list[str],
    save_path: Path,
) -> None:
    try:
        from bokeh.io import output_file, save
        from bokeh.models import CategoricalColorMapper, ColumnDataSource, HoverTool
        from bokeh.plotting import figure
    except ImportError as exc:
        print(f"Skipping interactive topography plot because Bokeh is unavailable: {exc}")
        return

    palette = palette_for_layers(layer_labels)
    hover_df = df.copy()
    hover_df["data_raster_png"] = build_hover_raster_previews(
        all_data=all_data,
        unit_ids_1based=hover_df["unit_id"].to_numpy(dtype=int),
        peak_indices=peak_indices,
    )
    hover_df["model_raster_png"] = build_model_hover_raster_previews(
        output=output,
        unit_ids_1based=hover_df["unit_id"].to_numpy(dtype=int),
        best_batch_indices_1based=hover_df["best_batch"].to_numpy(dtype=int),
    )

    source = ColumnDataSource(hover_df)
    color_mapper = CategoricalColorMapper(
        factors=layer_labels,
        palette=[palette[layer] for layer in layer_labels],
    )

    rgba_image, x0, y0, dw, dh = build_layer_topography_rgba(df, layer_labels)

    hover_html = """
    <div style="width: 690px;">
      <div><strong>Unit:</strong> @unit_id</div>
      <div><strong>Layer:</strong> @layer</div>
      <div><strong>Best batch:</strong> @best_batch</div>
      <div style="display: flex; gap: 12px; margin-top: 8px;">
        <div>
          <div style="font-weight: 600; margin-bottom: 4px;">Data raster</div>
          <img
            src="@data_raster_png"
            alt="Data raster preview"
            style="width: 320px; border: 1px solid #d0d0d0;"
          />
        </div>
        <div>
          <div style="font-weight: 600; margin-bottom: 4px;">Best sim raster</div>
          <img
            src="@model_raster_png"
            alt="Model raster preview"
            style="width: 320px; border: 1px solid #d0d0d0;"
          />
        </div>
      </div>
    </div>
    """

    plot = figure(
        title="UMAP projection with layer topography and data/model raster hover",
        width=950,
        height=760,
        tools="pan,wheel_zoom,box_zoom,reset,save,hover",
        active_scroll="wheel_zoom",
    )
    plot.image_rgba(
        image=[rgba_image],
        x=x0,
        y=y0,
        dw=dw,
        dh=dh,
        level="underlay",
    )
    points = plot.circle(
        x="umap_1",
        y="umap_2",
        source=source,
        size=9,
        alpha=0.95,
        line_color="white",
        line_width=0.6,
        fill_color={"field": "layer", "transform": color_mapper},
        legend_field="layer",
    )
    plot.xaxis.axis_label = "UMAP 1"
    plot.yaxis.axis_label = "UMAP 2"
    plot.legend.location = "top_right"
    plot.legend.click_policy = "hide"

    hover = plot.select_one(HoverTool)
    hover.tooltips = hover_html
    hover.point_policy = "follow_mouse"
    hover.renderers = [points]

    output_file(save_path, title="UMAP layer hover topography")
    save(plot)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading unit metadata from: {ALL_UNITS_PATH}")
    all_data = load_all_data(ALL_UNITS_PATH)

    print(f"Loading output, losses, and parameters from: {OUTPUT_MAT_PATH}")
    output, losses, params = load_output_bundle(OUTPUT_MAT_PATH)
    best_batch_indices_1based = best_batch_indices_from_final_epoch(losses)

    analysis_df, param_cols, layer_labels = build_analysis_dataframe(
        all_data=all_data,
        best_batch_indices_1based=best_batch_indices_1based,
        params=params,
    )
    analysis_df = compute_umap_embedding(analysis_df, param_cols)

    csv_path = OUTPUT_DIR / "best_fit_params_with_layers_and_umap.csv"
    analysis_df.to_csv(csv_path, index=False)
    print(f"Saved analysis table to: {csv_path}")

    pairplot_path = OUTPUT_DIR / "best_fit_param_pairplot_by_layer.png"
    save_parameter_pairplot(
        df=analysis_df,
        param_cols=param_cols,
        layer_labels=layer_labels,
        save_path=pairplot_path,
    )
    print(f"Saved seaborn pairplot to: {pairplot_path}")

    embedding_path = OUTPUT_DIR / "umap_embedding_by_layer_seaborn.png"
    save_seaborn_embedding_plot(
        df=analysis_df,
        layer_labels=layer_labels,
        save_path=embedding_path,
    )
    print(f"Saved seaborn UMAP plot to: {embedding_path}")

    if MAKE_BOKEH_HOVER:
        print("Computing peak directions and building hover raster previews...")
        peak_indices = compute_peak_direction_indices(all_data)
        hover_path = OUTPUT_DIR / "umap_embedding_by_layer_hover.html"
        save_bokeh_hover_embedding(
            df=analysis_df,
            all_data=all_data,
            peak_indices=peak_indices,
            layer_labels=layer_labels,
            save_path=hover_path,
        )
        print(f"Saved interactive Bokeh plot to: {hover_path}")

        hover_topography_path = OUTPUT_DIR / "umap_embedding_by_layer_hover_with_model_and_topography.html"
        save_bokeh_hover_embedding_with_model_and_topography(
            df=analysis_df,
            all_data=all_data,
            output=output,
            peak_indices=peak_indices,
            layer_labels=layer_labels,
            save_path=hover_topography_path,
        )
        print(f"Saved interactive Bokeh topography plot to: {hover_topography_path}")


if __name__ == "__main__":
    main()
