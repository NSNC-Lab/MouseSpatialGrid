from __future__ import annotations

import numpy as np
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.io.matlab import mat_struct

sns.set(style='white', context='notebook', rc={'figure.figsize':(14,10)})

OUTPUT_MAT_PATH = Path(
    r"C:\Users\ipboy\Documents\GitHub\ModelingEffort\Multi-Channel\Plotting\OliverDataPlotting\output_compressed_Eprop_All_cells_20260402_090929.mat"
)

def load_output_bundle(output_mat_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mat = loadmat(output_mat_path, variable_names=["output", "losses", "params"])

    if "output" not in mat or "losses" not in mat or "params" not in mat:
        keys = [k for k in mat.keys() if not k.startswith("__")]
        raise KeyError(
            f"Expected all of 'output', 'losses', and 'params' in {output_mat_path}. Available keys: {keys}"
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

    # Metric index 1 is the second saved loss term: the binned SSE / PSTH loss.
    final_epoch_binned_sse = losses[-1, 1, :, :]
    return np.argmin(final_epoch_binned_sse, axis=-1) + 1


def main() -> None:

    output, losses, params = load_output_bundle(OUTPUT_MAT_PATH)
    best_batch_indices_1based = best_batch_indices_from_final_epoch(losses)

    params_of_interest = []
    
    for k in range(220):
        params_of_interest.append(params[49,:,k,best_batch_indices_1based[k] - 1].tolist())  # Convert to 0-based indexing

    #penguins = pd.read_csv("https://raw.githubusercontent.com/allisonhorst/palmerpenguins/c19a904462482430170bfe2c718775ddb7dbb885/inst/extdata/penguins.csv")
    #penguins.head()

    #penguins = penguins.dropna()
    #penguins.species.value_counts()

    #sns.pairplot(penguins.drop("year", axis=1), hue='species')

    from umap import UMAP

    reducer = UMAP()

    #penguin_data = penguins[
    #    [
    #        "bill_length_mm",
    #        "bill_depth_mm",
    #        "flipper_length_mm",
    #        "body_mass_g",
    #    ]
    #].values



    scaled_params = StandardScaler().fit_transform(params_of_interest)

    embedding = reducer.fit_transform(scaled_params)
    embedding.shape

    plt.scatter(
        embedding[:, 0],
        embedding[:, 1],
        s=50)
    for i, (x_coord, y_coord) in enumerate(embedding):
        plt.text(x_coord, y_coord, str(i), fontsize=9, alpha=0.7)
        #c=[sns.color_palette()[x] for x in penguins.species.map({"Adelie":0, "Chinstrap":1, "Gentoo":2})])
    plt.gca().set_aspect('equal', 'datalim')
    plt.title('UMAP projection of the parameter space', fontsize=24)
    plt.show()


if __name__ == "__main__":
    main()