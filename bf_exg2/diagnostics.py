import csv
import os

import numpy as np
import matplotlib.pyplot as plt

from bf_exg2.exg2 import MyTrainer
from bf_exg2.default import DIRS

from bayesflow.diagnostics import plot_recovery


class Diagnostics:
    def __init__(self, trainer: MyTrainer, name=None) -> None:
        self.name = trainer.name if name is None else name
        self.generative_model = trainer.generative_model
        self.amortizer = trainer.amortizer

    def _validation_sim(self, batch_size, trials, **kwargs):
        prior_out = self.generative_model.prior(batch_size)
        n_stop = np.ceil(
            trials * self.generative_model.simulator.staircase["prop_ss"]
        ).astype(int)
        n_go = trials - n_stop
        context_batch = (
            trials,
            n_stop,
            n_go,
            self.generative_model.simulator.staircase["SSD_start"],
            self.generative_model.simulator.staircase["SSD_step"],
            self.generative_model.simulator.staircase["SSD_ceil"],
        )
        sim_data = self.generative_model.simulator.simulator(
            prior_out["prior_draws"], context_batch
        )
        out_dict = {
            "parameters": prior_out["prior_draws"],
            "summary_conditions": sim_data,
        }
        return out_dict

    def validation_sample(
        self,
        n_val_batch,
        n_samples,
        size=None,
        _return=False,
        **kwargs,
    ):
        if isinstance(self.generative_model.simulator.n_trials, int):
            validation_sims = self.generative_model.configurator(
                self.generative_model(batch_size=n_val_batch)
            )
            post_samples = self.amortizer.sample(validation_sims, n_samples=n_samples)
        elif isinstance(size, int):
            validation_sims = self._validation_sim(batch_size=n_val_batch, trials=size)
            print(validation_sims)
            post_samples = self.amortizer.sample(validation_sims, n_samples=n_samples)
        elif isinstance(size, dict):
            validation_sims = {}
            post_samples = {}
            for key, val in size.items():
                validation_sims[key] = self._validation_sim(
                    batch_size=n_val_batch, trials=val
                )
                post_samples[key] = self.amortizer.sample(
                    validation_sims[key], n_samples=n_samples
                )
        if _return:
            return validation_sims, post_samples
        else:
            self.validation_sims = validation_sims
            self.post_samples = post_samples

    def plot_posterior_recovery(self, standardized=False, **kwargs):
        if "validation_sims" in kwargs:
            validation_sims = kwargs["validation_sims"]
            post_samples = kwargs["post_samples"]
        else:
            validation_sims = self.validation_sims
            post_samples = self.post_samples
        if isinstance(post_samples, np.ndarray):
            if not standardized:
                post_samps = (
                    post_samples * self.generative_model.prior_sd
                    + self.generative_model.prior_mu
                )
                val_pars = (
                    validation_sims["parameters"] * self.generative_model.prior_sd
                    + self.generative_model.prior_mu
                )
            else:
                post_samps = post_samples
                val_pars = validation_sims["parameters"]
            plot_recovery(
                post_samps,
                val_pars,
                param_names=self.generative_model.prior.param_names,
            )
            if "save" in kwargs and kwargs["save"]:
                plt.savefig(f"{DIRS['figures']}/{self.name}.png")
        elif isinstance(post_samples, dict):
            for key in post_samples:
                if not standardized:
                    post_samps = (
                        post_samples[key] * self.generative_model.prior_sd
                        + self.generative_model.prior_mu
                    )
                    val_pars = (
                        validation_sims[key]["parameters"]
                        * self.generative_model.prior_sd
                        + self.generative_model.prior_mu
                    )
                else:
                    post_samps = post_samples[key]
                    val_pars = validation_sims[key]["parameters"]
                plot_recovery(
                    post_samps,
                    val_pars,
                    param_names=self.generative_model.prior.param_names,
                )
                if "save" in kwargs and kwargs["save"]:
                    plt.savefig(f"{DIRS['figures']}/{self.name}_{key}.png")


class Inference:
    def __init__(self, trainer: MyTrainer, name=None) -> None:
        self.name = trainer.name if name is None else name
        self.generative_model = trainer.generative_model
        self.amortizer = trainer.amortizer

    def sample(self, n_data_sets, size=None, save=True, **kwargs):
        if isinstance(self.generative_model.simulator.n_trials, int):
            data = self.generative_model.configurator(
                self.generative_model(batch_size=n_data_sets)
            )
        if isinstance(self.generative_model.simulator.n_trials, tuple):
            if size is None:
                params = np.zeros(
                    (n_data_sets, len(self.generative_model.prior.param_names))
                )
                sim_data = np.zeros(
                    (n_data_sets, self.generative_model.simulator.n_trials[1], 2)
                )
                for i in range(n_data_sets):
                    d = self.generative_model.configurator(
                        self.generative_model(batch_size=1)
                    )
                    params[i, :] = d["parameters"][0]
                    sim_data[i, :, :] = d["summary_conditions"][0]
                data = {"parameters": params, "summary_conditions": sim_data}
            elif isinstance(size, int):
                prior_out = self.generative_model.prior(n_data_sets)
                n_stop = np.ceil(
                    size * self.generative_model.simulator.staircase["prop_ss"]
                ).astype(int)
                n_go = size - n_stop
                context_batch = (
                    size,
                    n_stop,
                    n_go,
                    self.generative_model.simulator.staircase["SSD_start"],
                    self.generative_model.simulator.staircase["SSD_step"],
                    self.generative_model.simulator.staircase["SSD_ceil"],
                )
                sim_data = self.generative_model.simulator.simulator(
                    prior_out["prior_draws"], context_batch
                )
                data = {
                    "parameters": prior_out["prior_draws"],
                    "summary_conditions": sim_data,
                }
        if save:
            self.parameters = data["parameters"]
            self.data = data["summary_conditions"]
        return data

    def save_for_dmc(self, name=None):
        if name is None:
            name = self.name
        os.makedirs(f"{DIRS['validations_data']}/{name}", exist_ok=True)
        for i in range(self.parameters.shape[0]):
            with open(
                f"{DIRS['validations_data']}/{name}/{name}_{i}.csv",
                "wt",
                newline="",
                encoding="utf-8",
            ) as f:
                writer = csv.writer(f)
                writer.writerows(self.data[i, :, :])
        with open(
            f"{DIRS['validations_data']}/{name}_params.csv",
            "wt",
            newline="",
            encoding="utf-8",
        ) as f:
            writer = csv.writer(f)
            writer.writerows(
                self.parameters * self.generative_model.prior_sd
                + self.generative_model.prior_mu
            )

    def sample_posterior(self, n_samples, _return=False, **kwargs):
        if "data" in kwargs:
            data = kwargs["data"]
        else:
            data = self.data
        if "parameters" in kwargs:
            parameters = kwargs["parameters"]
        else:
            parameters = self.parameters

        use_data = {
            "parameters": parameters,
            "summary_conditions": data,
        }
        posterior_samples = self.amortizer.sample(use_data, n_samples=n_samples)
        self.posterior_samples = posterior_samples
        if _return:
            return posterior_samples
