from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


def plot_noise_raster(
    noise: np.ndarray,
    *,
    cell: int,
    batch: int,
    channel: int = 1,
    dt_ms: float = 1.0,
    one_indexed: bool = True,
    title: str | None = None,
) -> None:
    """Plot a raster for a selected (cell, batch) slice of `noise`.

    Expected shape (as used in run_main_integrated_all_cells.py):
        (num_cells, batch, trials, time, chans)

    Args:
        noise: spike tensor, typically uint8/bool.
        cell: which cell to view (1-indexed by default to match MATLAB).
        batch: which batch item to view (1-indexed by default to match MATLAB).
        channel: which channel to view (1-indexed by default).
        dt_ms: timestep in ms for x-axis scaling.
        one_indexed: interpret indices as MATLAB-style if True.
        title: optional custom title.
    """

    if noise.ndim != 5:
        raise ValueError(f"Expected noise.ndim == 5, got shape {noise.shape}")

    cell_idx = cell - 1 if one_indexed else cell
    batch_idx = batch - 1 if one_indexed else batch
    chan_idx = channel - 1 if one_indexed else channel

    if not (0 <= cell_idx < noise.shape[0]):
        raise IndexError(f"cell out of range: {cell} for noise.shape[0]={noise.shape[0]}")
    if not (0 <= batch_idx < noise.shape[1]):
        raise IndexError(f"batch out of range: {batch} for noise.shape[1]={noise.shape[1]}")
    if not (0 <= chan_idx < noise.shape[4]):
        raise IndexError(f"channel out of range: {channel} for noise.shape[4]={noise.shape[4]}")

    # trials x time
    trial_time = noise[cell_idx, batch_idx, :, :, chan_idx]
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
        plt.title(f"Noise raster: cell={cell}, batch={batch}, channel={channel}")
    else:
        plt.title(title)
    plt.tight_layout()
    plt.show()
