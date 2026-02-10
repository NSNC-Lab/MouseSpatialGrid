from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


def plot_data_raster(
    data: np.ndarray,
    *,
    cell: int,
    channel: int = 1,
    dt_ms: float = 1.0,
    one_indexed: bool = True,
    title: str | None = None,
) -> None:
    """Plot a raster for a selected cell from `data`.

    Expected shapes:
        - (cells, trials, time, chans)
        - (cells, trials, time, 1)

    This matches your described `data` shape: (220, 10, 29801, 1).

    Args:
        data: spike tensor, typically uint8/bool.
        cell: which cell to view (1-indexed by default to match MATLAB).
        channel: which channel to view (1-indexed by default).
        dt_ms: timestep in ms for x-axis scaling.
        one_indexed: interpret indices as MATLAB-style if True.
        title: optional custom title.
    """

    if data.ndim != 4:
        raise ValueError(f"Expected data.ndim == 4, got shape {data.shape}")

    cell_idx = cell - 1 if one_indexed else cell
    chan_idx = channel - 1 if one_indexed else channel

    if not (0 <= cell_idx < data.shape[0]):
        raise IndexError(f"cell out of range: {cell} for data.shape[0]={data.shape[0]}")
    if not (0 <= chan_idx < data.shape[3]):
        raise IndexError(
            f"channel out of range: {channel} for data.shape[3]={data.shape[3]}"
        )

    # trials x time
    trial_time = data[cell_idx, :, :, chan_idx]
    if trial_time.ndim != 2:
        trial_time = np.squeeze(trial_time)
    if trial_time.ndim != 2:
        raise ValueError(
            f"Expected selected slice to be 2D (trials,time), got {trial_time.shape}"
        )

    spike_times_by_trial = []
    for tr in range(trial_time.shape[0]):
        spike_inds = np.flatnonzero(trial_time[tr] > 0)
        spike_times_by_trial.append(spike_inds.astype(float) * float(dt_ms))

    plt.figure(figsize=(10, 4))
    plt.eventplot(spike_times_by_trial, colors="k", linewidths=0.8)
    plt.xlabel("Time (ms)")
    plt.ylabel("Trial")
    if title is None:
        plt.title(f"Data raster: cell={cell}, channel={channel}")
    else:
        plt.title(title)
    plt.tight_layout()
    plt.show()
