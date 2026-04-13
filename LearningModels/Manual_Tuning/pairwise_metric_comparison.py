from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_HISTORY_PATH = Path(
    r"C:\Users\ipboy\Documents\GitHub\MouseSpatialGrid\LearningModels\Manual_Tuning\post_fit_reports\unit_133\post_fit_unit_history.csv"
)
DEFAULT_TARGET_METRIC = "mds_manual_fit_distance_to_data"


def main() -> None:
    args = parse_args()
    history_path = args.history.resolve()
    if not history_path.exists():
        raise FileNotFoundError(f"Could not find history CSV: {history_path}")

    output_dir = (args.output_dir or history_path.parent).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = read_csv_rows(history_path)
    if not rows:
        raise ValueError(f"No rows found in {history_path}")

    numeric_columns = find_numeric_columns(rows, min_count=args.min_count)
    if args.target_metric not in numeric_columns:
        raise ValueError(
            f"Target metric '{args.target_metric}' is not available as a numeric column. "
            f"Numeric columns found: {', '.join(numeric_columns)}"
        )

    data = build_numeric_arrays(rows, numeric_columns)
    target_rows = compute_target_correlations(
        data=data,
        target_metric=args.target_metric,
        columns=numeric_columns,
        min_count=args.min_count,
    )
    matrix_rows, long_rows = compute_pairwise_correlations(
        data=data,
        columns=numeric_columns,
        min_count=args.min_count,
    )

    target_csv_path = output_dir / f"{history_path.stem}_{sanitize(args.target_metric)}_pearson.csv"
    matrix_csv_path = output_dir / f"{history_path.stem}_pairwise_pearson_matrix.csv"
    long_csv_path = output_dir / f"{history_path.stem}_pairwise_pearson_long.csv"
    summary_txt_path = output_dir / f"{history_path.stem}_{sanitize(args.target_metric)}_summary.txt"
    heatmap_svg_path = output_dir / f"{history_path.stem}_pairwise_pearson_heatmap.svg"

    target_csv_path = write_dict_rows_csv(target_csv_path, target_rows)
    matrix_csv_path = write_matrix_csv(matrix_csv_path, numeric_columns, matrix_rows)
    long_csv_path = write_dict_rows_csv(long_csv_path, long_rows)
    summary_txt_path = write_summary_txt(
        summary_txt_path,
        history_path=history_path,
        target_metric=args.target_metric,
        target_rows=target_rows,
        pairwise_rows=long_rows,
    )
    heatmap_svg_path = save_pairwise_heatmap(
        path=heatmap_svg_path,
        columns=numeric_columns,
        matrix_rows=matrix_rows,
        title=f"Pairwise Pearson Correlations\n{history_path.stem}",
    )

    print(f"Loaded history rows: {len(rows)}")
    print(f"Numeric columns analyzed: {len(numeric_columns)}")
    print(f"Saved target-vs-all correlations: {target_csv_path}")
    print(f"Saved pairwise matrix: {matrix_csv_path}")
    print(f"Saved pairwise long table: {long_csv_path}")
    print(f"Saved summary text: {summary_txt_path}")
    print(f"Saved pairwise heatmap: {heatmap_svg_path}")

    print("\nTop 10 |r| values against the target metric:")
    for row in target_rows[:10]:
        print(
            f"  {row['metric']}: r={row['pearson_r']}, n={row['pair_count']}, "
            f"abs_r={row['abs_pearson_r']}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute correlations from a manual-tuning post-fit history CSV. "
            "Outputs target-vs-all Pearson correlations plus a full pairwise Pearson matrix."
        )
    )
    parser.add_argument(
        "--history",
        type=Path,
        default=DEFAULT_HISTORY_PATH,
        help="Path to a post_fit_unit_history.csv file.",
    )
    parser.add_argument(
        "--target-metric",
        default=DEFAULT_TARGET_METRIC,
        help="Metric to compare against all other numeric metrics.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for correlation outputs. Defaults to the history file's folder.",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=3,
        help="Minimum number of paired finite observations required to report a correlation.",
    )
    return parser.parse_args()


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def find_numeric_columns(rows: Sequence[Dict[str, str]], min_count: int) -> List[str]:
    if not rows:
        return []

    numeric_columns: List[str] = []
    for column in rows[0].keys():
        finite_values = []
        for row in rows:
            value = parse_float(row.get(column, ""))
            if value is not None and np.isfinite(value):
                finite_values.append(value)
        if len(finite_values) >= min_count:
            numeric_columns.append(column)
    return numeric_columns


def build_numeric_arrays(rows: Sequence[Dict[str, str]], columns: Sequence[str]) -> Dict[str, np.ndarray]:
    arrays: Dict[str, np.ndarray] = {}
    for column in columns:
        values = [parse_float(row.get(column, "")) for row in rows]
        arrays[column] = np.array(
            [value if value is not None else np.nan for value in values],
            dtype=np.float64,
        )
    return arrays


def compute_target_correlations(
    *,
    data: Dict[str, np.ndarray],
    target_metric: str,
    columns: Sequence[str],
    min_count: int,
) -> List[Dict[str, object]]:
    target = data[target_metric]
    rows: List[Dict[str, object]] = []

    for column in columns:
        if column == target_metric:
            continue
        r_value, pair_count = pearson_correlation(target, data[column], min_count=min_count)
        rows.append(
            {
                "target_metric": target_metric,
                "metric": column,
                "pearson_r": r_value,
                "abs_pearson_r": abs(r_value) if np.isfinite(r_value) else np.nan,
                "pair_count": pair_count,
            }
        )

    rows.sort(
        key=lambda row: (
            not np.isfinite(row["abs_pearson_r"]),
            -(row["abs_pearson_r"] if np.isfinite(row["abs_pearson_r"]) else 0.0),
            row["metric"],
        )
    )
    return rows


def compute_pairwise_correlations(
    *,
    data: Dict[str, np.ndarray],
    columns: Sequence[str],
    min_count: int,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    matrix_rows: List[Dict[str, object]] = []
    long_rows: List[Dict[str, object]] = []

    for row_index, row_metric in enumerate(columns):
        row: Dict[str, object] = {"metric": row_metric}
        for col_index, col_metric in enumerate(columns):
            r_value, pair_count = pearson_correlation(data[row_metric], data[col_metric], min_count=min_count)
            row[col_metric] = r_value
            if row_index < col_index:
                long_rows.append(
                    {
                        "metric_a": row_metric,
                        "metric_b": col_metric,
                        "pearson_r": r_value,
                        "abs_pearson_r": abs(r_value) if np.isfinite(r_value) else np.nan,
                        "pair_count": pair_count,
                    }
                )
        matrix_rows.append(row)

    long_rows.sort(
        key=lambda row: (
            not np.isfinite(row["abs_pearson_r"]),
            -(row["abs_pearson_r"] if np.isfinite(row["abs_pearson_r"]) else 0.0),
            row["metric_a"],
            row["metric_b"],
        )
    )
    return matrix_rows, long_rows


def pearson_correlation(x: np.ndarray, y: np.ndarray, *, min_count: int) -> Tuple[float, int]:
    mask = np.isfinite(x) & np.isfinite(y)
    pair_count = int(np.sum(mask))
    if pair_count < min_count:
        return np.nan, pair_count

    x_valid = x[mask]
    y_valid = y[mask]
    if np.allclose(x_valid, x_valid[0]) or np.allclose(y_valid, y_valid[0]):
        return np.nan, pair_count

    return float(np.corrcoef(x_valid, y_valid)[0, 1]), pair_count


def parse_float(raw_value: str | None) -> float | None:
    if raw_value is None:
        return None
    text = str(raw_value).strip()
    if text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def write_matrix_csv(path: Path, columns: Sequence[str], rows: Iterable[Dict[str, object]]) -> Path:
    headers = ["metric", *columns]
    return write_csv_rows(path, headers, rows)


def write_dict_rows_csv(path: Path, rows: Iterable[Dict[str, object]]) -> Path:
    rows = list(rows)
    if not rows:
        return write_text_with_fallback(path, "")
    headers = list(rows[0].keys())
    return write_csv_rows(path, headers, rows)


def write_csv_rows(path: Path, headers: Sequence[str], rows: Iterable[Dict[str, object]]) -> Path:
    def _writer(output_path: Path) -> None:
        with output_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(headers), extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow({header: stringify_value(row.get(header, "")) for header in headers})

    return write_with_fallback(path, _writer)


def write_summary_txt(
    path: Path,
    *,
    history_path: Path,
    target_metric: str,
    target_rows: Sequence[Dict[str, object]],
    pairwise_rows: Sequence[Dict[str, object]],
) -> Path:
    top_target = target_rows[:10]
    top_pairwise = pairwise_rows[:15]
    lines = [
        f"History file: {history_path}",
        f"Target metric: {target_metric}",
        "",
        "Top 10 absolute Pearson correlations against the target metric:",
    ]
    for row in top_target:
        lines.append(
            f"  {row['metric']}: r={row['pearson_r']}, n={row['pair_count']}, abs_r={row['abs_pearson_r']}"
        )

    lines.append("")
    lines.append("Top 15 absolute Pearson correlations across all metric pairs:")
    for row in top_pairwise:
        lines.append(
            f"  {row['metric_a']} vs {row['metric_b']}: "
            f"r={row['pearson_r']}, n={row['pair_count']}, abs_r={row['abs_pearson_r']}"
        )

    return write_text_with_fallback(path, "\n".join(lines) + "\n")


def save_pairwise_heatmap(
    *,
    path: Path,
    columns: Sequence[str],
    matrix_rows: Sequence[Dict[str, object]],
    title: str,
) -> Path:
    matrix = matrix_rows_to_array(columns, matrix_rows)
    n_metrics = len(columns)
    if n_metrics == 0:
        return path

    fig_size = max(10.0, min(0.42 * n_metrics + 6.0, 24.0))
    label_fontsize = max(5, min(10, int(round(220 / max(n_metrics, 1)))))

    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad(color="#d9d9d9")

    image = ax.imshow(matrix, cmap=cmap, vmin=-1.0, vmax=1.0, aspect="equal")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Pearson r")

    ax.set_xticks(np.arange(n_metrics))
    ax.set_yticks(np.arange(n_metrics))
    ax.set_xticklabels(columns, rotation=90, fontsize=label_fontsize)
    ax.set_yticklabels(columns, fontsize=label_fontsize)
    ax.set_title(title)
    ax.set_xlabel("Metric")
    ax.set_ylabel("Metric")

    ax.set_xticks(np.arange(-0.5, n_metrics, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_metrics, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.4, alpha=0.7)
    ax.tick_params(which="minor", bottom=False, left=False)

    fig.tight_layout()
    saved_path = save_figure_with_fallback(fig, path, format="svg", bbox_inches="tight")
    plt.close(fig)
    return saved_path


def matrix_rows_to_array(columns: Sequence[str], matrix_rows: Sequence[Dict[str, object]]) -> np.ndarray:
    matrix = np.full((len(columns), len(columns)), np.nan, dtype=np.float64)
    for row_index, row in enumerate(matrix_rows):
        for col_index, column in enumerate(columns):
            value = parse_float(row.get(column))
            if value is not None:
                matrix[row_index, col_index] = value
    return matrix


def stringify_value(value: object) -> str:
    if isinstance(value, (float, np.floating)):
        if np.isnan(value):
            return "nan"
        return f"{float(value):.12g}"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    return str(value)


def sanitize(name: str) -> str:
    chars = []
    for char in str(name):
        if char.isalnum():
            chars.append(char)
        else:
            chars.append("_")
    compact = "".join(chars).strip("_")
    while "__" in compact:
        compact = compact.replace("__", "_")
    return compact


def write_text_with_fallback(path: Path, text: str) -> Path:
    def _writer(output_path: Path) -> None:
        output_path.write_text(text, encoding="utf-8")

    return write_with_fallback(path, _writer)


def save_figure_with_fallback(fig, path: Path, **savefig_kwargs) -> Path:
    def _writer(output_path: Path) -> None:
        fig.savefig(output_path, **savefig_kwargs)

    return write_with_fallback(path, _writer)


def write_with_fallback(path: Path, writer) -> Path:
    try:
        writer(path)
        return path
    except PermissionError:
        fallback = timestamped_fallback_path(path)
        print(f"Output file is locked: {path}")
        print(f"Saving instead to: {fallback}")
        writer(fallback)
        return fallback


def timestamped_fallback_path(path: Path) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return path.with_name(f"{path.stem}_{timestamp}{path.suffix}")


if __name__ == "__main__":
    main()
