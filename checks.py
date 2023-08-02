from bf_exg2.exg2 import MyPrior, MySimulator, MyGenerativeModel, MyTrainer
from bf_exg2.diagnostics import Diagnostics
from bf_exg2.default import PARAM_NAMES, DIRS, default_prior
from bf_exg2.transform import boxcox_inv as untrans

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import pickle as pkl

SIZE = 1000

NAME = "checks"

recov_training = {
    "num_epochs": 128,
    "iter_per_epoch": 512,
    "batch_size": 64,
}


def plot_prior(prior_samples, post_samples):
    f, axarr = plt.subplots(3, 2, figsize=(18, 8))
    for i, ax in enumerate(axarr.flat):
        sns.kdeplot(
            prior_samples[:, i],
            ax=ax,
            fill=True,
            color="black",
            alpha=0.5,
            label="Prior distribution",
        )
        sns.kdeplot(
            post_samples[:, i],
            ax=ax,
            fill=True,
            color="maroon",
            alpha=0.5,
            label="Posterior distribution",
        )
        # ax.set_xlim((0, 4))
        ax.set_xlabel("Parameter values")
        ax.set_title(PARAM_NAMES[i])
        if i == 1:
            ax.legend()
        sns.despine()
        f.tight_layout()

    plt.savefig(f"{DIRS['figures']}/prior_recovery.png", dpi=300)


def plot_rt_distribution(data, bins=50):
    unique_ssd = np.unique(data[data[:, 1] != -1.0, 1])
    multi_rt = [np.setdiff1d(data[data[:, 1] == ssd, 0], -1.0) for ssd in unique_ssd]
    cols = mpl.colormaps["viridis"](np.linspace(0, 1, len(unique_ssd)))

    ppf_plot, ax = plt.subplots(1, 2, figsize=(12, 4))
    ppf_plot.suptitle("Signal Response Time Distribution")
    ax[0].hist(  # with elapsed responses
        data[data[:, 1] == -1.0, 0],
        # ax=ax[0],
        fill=True,
        color=cols[0],
        density=True,
        stacked=True,
        bins=bins,
        label="Go RTs",
    )
    ax[1].hist(  # with elapsed responses
        np.setdiff1d(data[data[:, 1] != -1.0, 0], -1.0),
        # ax=ax[1],
        fill=True,
        color=cols[0],
        density=True,
        bins=bins,
        label="All Stop RTs",
        alpha=0.5,
    )
    for i, ssd in enumerate(unique_ssd):
        if i > 2 and i < len(unique_ssd) - 2:
            z, temp_bins = np.histogram(multi_rt[i], bins=bins)
            bincenters = 0.5 * (temp_bins[1:] + temp_bins[:-1])
            ax[1].plot(
                bincenters,
                z / np.sum(z) * 40,
                color=cols[-i],
                label=f"{round(ssd * 1000)}ms",
            )
    ax[1].legend(title="SSD")
    ax[0].set_xlim((0, 1.0))
    ax[1].set_xlim((0, 1.0))
    ax[0].set_ylim((0, 5.5))
    ax[1].set_ylim((0, 7.0))
    ax[0].set_xlabel("Go RT (s)")
    ax[1].set_xlabel("Stop RT (s)")
    ax[0].set_yticklabels([])
    ax[1].set_yticklabels([])

    plt.savefig(f"{DIRS['figures']}/ppf_rt-dist.png", dpi=900)


def main():
    prior = MyPrior("mvn", lhs_sampling=True)
    mvn_prior = prior.plot_prior2d()
    plt.savefig(f"{DIRS['figures']}/mvn_prior.png", dpi=300)
    prio2 = MyPrior("uvn", lhs_sampling=True)
    uvn_prior = prio2.plot_prior2d()
    plt.savefig(f"{DIRS['figures']}/uvn_prior.png", dpi=300)
    simulator = MySimulator(n_trials=SIZE, uninformative=True)
    model = MyGenerativeModel(prior, simulator)
    trainer = MyTrainer(model, name=NAME)
    trainer.train_online(training_settings=recov_training)
    # pkl.dump(loss, open(f"{NAME}.pkl", "wb"))
    diagnostics = Diagnostics(trainer)
    _, post_samples = diagnostics.validation_sample(1, 10000, SIZE, _return=True)
    prior_samples = prior(10000)["prior_draws"]
    plot_prior(prior_samples, post_samples * model.prior_sd + model.prior_mu)

    ppf_simulator = MySimulator(n_trials=SIZE, uninformative=False)
    prior_vals = default_prior()
    ppf_data = ppf_simulator(
        untrans(np.asarray([prior_vals["mu"]]), prior_vals["trans"])
    )["sim_data"][0]

    plot_rt_distribution(ppf_data)

    # sns.despine()


if __name__ == "__main__":
    main()
