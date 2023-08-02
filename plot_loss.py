import os
import pickle as pkl

import matplotlib.pyplot as plt
import numpy as np

from bf_exg2.default import DIRS, FULL_TRAINING_SETTINGS

MODELS = {
    "mvn-lhs": "multivariate normal, LHS-sampled",
    "mvn-rand": "multivariate normal, random sampling",
    "uvn-lhs": "univariate normal, LHS-sampled",
    "uvn-rand": "univariate normal, random sampling",
    "uvu-lhs-trans": "univariate uniform, LHS-sampled from transformed space",
    "uvu-rand-trans": "univariate uniform, random sampling from transformed space",
    "uvu-lhs-untr": "uniform, LHS-sampled from untransformed space",
    "uvu-rand-untr": "uniform, random sampling from untransformed space",
    "mvn-lhs_10k": "MVN fixed at 10k observations",
    "uvn-lhs_10k": "UVN fixed at 10k observations",
}


def calculate_mean_std(array, n):
    array_length = len(array)
    num_chunks = array_length // n  # Number of chunks of size n
    remainder = array_length % n  # Remaining elements

    # Reshape the array to divide it into chunks of size n
    reshaped_array = array[: array_length - remainder].reshape(num_chunks, n)

    # Calculate the mean and standard deviation for each chunk
    chunk_means = np.mean(reshaped_array, axis=1)
    chunk_stds = np.std(reshaped_array, axis=1)

    # If there are remaining elements, calculate the mean and standard deviation for them
    if remainder > 0:
        remaining_elements = array[-remainder:]
        remaining_mean = np.mean(remaining_elements)
        remaining_std = np.std(remaining_elements)
        chunk_means = np.append(chunk_means, remaining_mean)
        chunk_stds = np.append(chunk_stds, remaining_std)

    return np.stack((chunk_means, chunk_stds), axis=1)


def main():
    loss_pkl = os.listdir("loss")
    loss_dict = {}
    for file in loss_pkl:
        loss_dict[file[:-4]] = pkl.load(open(f"loss/{file}", "rb"))
    plt.figure(figsize=(10, 5))
    for key, value in loss_dict.items():
        # if key in ["uvu-lhs-trans", "uvu-lhs-untr", "uvu-rand-trans", "uvu-rand-untr"]:
        if key in [
            "mvn-lhs",
            "mvn-rand",
            "uvn-lhs",
            "uvn-rand",
            "uvu-rand-untr",
            "uvu-lhs-untr",
        ]:
            # if key in ["mvn-lhs_10k", "uvn-lhs_10k"]:
            mean_std = calculate_mean_std(
                np.asarray(value.iloc[:, 0]),
                FULL_TRAINING_SETTINGS["iter_per_epoch"] // 4,
            )
            plt.errorbar(
                np.arange(len(mean_std)),
                mean_std[:, 0],
                yerr=mean_std[
                    :, 1
                ],  # / np.sqrt(FULL_TRAINING_SETTINGS["iter_per_epoch"]),
                label=MODELS[key],
                elinewidth=0.05,
            )
    plt.ylim(-4, 3)
    plt.title("Average Loss per Epoch")
    plt.xlabel(f"Epoch ({FULL_TRAINING_SETTINGS['iter_per_epoch'] // 4} it./ep.)")
    plt.ylabel("Default Loss (w/ SD per ep.)")
    plt.legend(loc="upper right")
    plt.savefig(f"{DIRS['figures']}/loss_normal_per_sampling.png", dpi=900)
