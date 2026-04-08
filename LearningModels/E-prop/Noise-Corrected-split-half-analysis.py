import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Point this at your local NEMS correlation.py
corr_path = Path(r"C:\Users\ipboy\Documents\New project\NEMS_repo\nems\metrics\correlation.py")

spec = importlib.util.spec_from_file_location("nems_corr", corr_path)
nems_corr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(nems_corr)

correlation = nems_corr.correlation
noise_corrected_r = nems_corr.noise_corrected_r

rng = np.random.default_rng(0)

# Fake example data
# single_trial_responses shape: (trials, time, neurons)
# prediction shape: (time, neurons)
n_trials = 20
n_time = 1000
n_neurons = 12

true_signal = rng.normal(size=(n_time, n_neurons))
single_trial_responses = true_signal[None, :, :] + 0.8 * rng.normal(
    size=(n_trials, n_time, n_neurons)
)
prediction = true_signal + 0.4 * rng.normal(size=(n_time, n_neurons))

# Raw correlation: model vs trial-averaged response
mean_response = single_trial_responses.mean(axis=0)
raw_r = correlation(prediction, mean_response)

# Noise-corrected correlation: model vs data
nc_r = noise_corrected_r(
    prediction,
    single_trial_responses,
    trial_axis=0,
    channel_axis=-1,
)

# Split-half reliability: data vs data
trial_order = rng.permutation(n_trials)
half1 = single_trial_responses[trial_order[: n_trials // 2]].mean(axis=0)
half2 = single_trial_responses[trial_order[n_trials // 2 :]].mean(axis=0)
split_half_r = correlation(half1, half2)

# Optional Spearman-Brown correction to estimate full-data reliability
split_half_sb = (2 * split_half_r) / (1 + split_half_r)

print("raw_r:", raw_r)
print("noise_corrected_r:", nc_r)
print("split_half_r:", split_half_r)
print("split_half_spearman_brown:", split_half_sb)

# Scatter: reliability vs model performance
plt.figure(figsize=(6, 6))
plt.scatter(split_half_sb, nc_r, alpha=0.8)
plt.xlabel("Data reliability (split-half, Spearman-Brown corrected)")
plt.ylabel("Model performance (noise-corrected r)")
plt.title("Model vs Data Reliability")
plt.axline((0, 0), slope=1, linestyle="--", color="gray")
plt.show()