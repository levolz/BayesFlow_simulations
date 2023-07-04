import numpy as np
import matplotlib.pyplot as plt

from bf_exg2 import sampling
from bf_exg2.model import simulate_batch, uninf
from bf_exg2.default import (
    PARAM_NAMES,
    DIRS,
    STAIRCASE,
    TRAINER_SETUP,
    TRAINING_SETTINGS,
    FULL_TRAINING_SETTINGS,
    default_prior,
)
from bf_exg2.transform import retransform

from bayesflow.simulation import Prior, ContextGenerator, Simulator, GenerativeModel
from bayesflow.networks import InvariantNetwork, InvertibleNetwork
from bayesflow.amortizers import AmortizedPosterior
from bayesflow.trainers import Trainer

from bayesflow.diagnostics import plot_losses


class MyPrior(Prior):
    def __init__(
        self,
        dist: str,
        params: dict = None,
        lhs_sampling: bool = True,
        trans: bool = True,
        **kwargs,
    ):
        self.param_names = (
            kwargs["param_names"] if "param_names" in kwargs else PARAM_NAMES
        )

        if params is None:
            set_uvu = True
            params = default_prior()
            trans_vec = params["trans"]
            bounds = (params["lower"], params["upper"])
        else:
            set_uvu = False

        if dist == "uvu":
            if set_uvu:
                if trans:
                    self.args = (
                        params["mu"] - 6 * params["sd"],
                        params["mu"] + 6 * params["sd"],
                    )
                else:
                    self.args = (
                        retransform(
                            np.expand_dims(params["lower"], axis=0), bounds, trans_vec
                        )[0],
                        retransform(
                            np.expand_dims(params["upper"], axis=0), bounds, trans_vec
                        )[0],
                    )
            else:
                self.args = (params["lower"], params["upper"])

            if lhs_sampling:
                self.dist = sampling.uvu_lhs
            else:
                self.dist = sampling.uvu
        elif dist == "uvn":
            self.args = (params["mu"], params["sd"])
            if lhs_sampling:
                self.dist = sampling.uvn_lhs
            else:
                self.dist = sampling.uvn
        elif dist == "mvn":
            self.args = (params["mu"], params["sigma"])
            if lhs_sampling:
                self.dist = sampling.mvn_lhs
            else:
                self.dist = sampling.mvn

        @staticmethod
        def batch_prior_fun(batch_size, args=self.args):
            if not trans:
                return self.dist(size=batch_size, *args)
            else:
                return retransform(self.dist(size=batch_size, *args), bounds, trans_vec)

        self.batched_prior = batch_prior_fun
        self.context_gen = None
        self.is_batched = True


class MySimulator(Simulator):
    staircase = STAIRCASE

    def __init__(self, n_trials, uninformative=False, **kwargs):
        if "SSD_start" in kwargs:
            self.staircase["SSD_start"] = kwargs["SSD_start"]
        if "SSD_step" in kwargs:
            self.staircase["SSD_step"] = kwargs["SSD_step"]
        if "SSD_ceil" in kwargs:
            self.staircase["SSD_ceil"] = kwargs["SSD_ceil"]
        if "prop_ss" in kwargs:
            self.staircase["prop_ss"] = kwargs["prop_ss"]

        self.n_trials = n_trials

        if uninformative:
            self.simulator = uninf
        else:
            self.simulator = simulate_batch
        self.is_batched = True

        @staticmethod
        def context_generator(
            n_trials=self.n_trials,
            staircase=self.staircase,
        ):
            n_trials = MySimulator.trials(n_trials)
            n_stop = np.ceil(n_trials * staircase["prop_ss"]).astype(int)
            n_go = n_trials - n_stop
            return (
                n_trials,
                n_stop,
                n_go,
                staircase["SSD_start"],
                staircase["SSD_step"],
                staircase["SSD_ceil"],
            )

        self.context_gen = ContextGenerator(non_batchable_context_fun=context_generator)

    @staticmethod
    def trials(n):
        if isinstance(n, int) or isinstance(n, float):
            out = n
        elif isinstance(n, tuple) and len(n) == 2:
            out = np.random.randint(*n)
        return out


class MyGenerativeModel(GenerativeModel):
    _N_STANDARD = 10000

    def __init__(self, prior, simulator, **kwargs):
        super().__init__(prior, simulator, **kwargs)

        self.model_name = kwargs["model_name"] if "model_name" in kwargs else "exg2"

        self.n_draws = kwargs["n_draws"] if "n_draws" in kwargs else self._N_STANDARD

        prior_mu, prior_sd = self.prior.estimate_means_and_stds(n_draws=self.n_draws)
        self.prior_mu = prior_mu
        self.prior_sd = prior_sd

        self.num_params = len(self.prior.param_names)

    def plot_pushforward(
        self,
        parameter_draws=None,
        funcs_list=None,
        funcs_labels=None,
        batch_size=1000,
        show_raw_sims=True,
    ):
        super().plot_pushforward(
            parameter_draws, funcs_list, funcs_labels, batch_size, show_raw_sims
        )

    def configurator(self, sim_dict, **kwargs):
        out = dict()

        data = sim_dict["sim_data"].astype(np.float32)
        out["summary_conditions"] = data

        params = sim_dict["prior_draws"].astype(np.float32)
        out["parameters"] = (params - self.prior_mu) / self.prior_sd

        return out


class MyTrainer(Trainer):
    trainer_setup = TRAINER_SETUP
    train_settings = TRAINING_SETTINGS
    checkpoint_path = DIRS["checkpoints"]
    figures_path = DIRS["figures"]

    def __init__(self, generative_model, name="anon", **kwargs):
        if isinstance(generative_model, MyGenerativeModel):
            if "summary_dim" in kwargs:
                self.trainer_setup["summary_dim"] = kwargs["summary_dim"]
            if "num_coupling_layers" in kwargs:
                self.trainer_setup["num_coupling_layers"] = kwargs[
                    "num_coupling_layers"
                ]
            if "checkpoint_path" in kwargs:
                self.checkpoint_path = kwargs["checkpoint_path"]
            if "figures_path" in kwargs:
                self.figures_path = kwargs["figures_path"]

            self.name = name
            self.generative_model = generative_model

            # if "checkpoints" in kwargs:
            #     checkpoints = kwargs["checkpoints"]
            # else:
            #     checkpoints = DIRS["checkpoints"]

            self.amortizer = AmortizedPosterior(
                InvertibleNetwork(
                    num_params=generative_model.num_params,
                    num_coupling_layers=self.trainer_setup["num_coupling_layers"],
                ),
                InvariantNetwork(summary_dim=self.trainer_setup["summary_dim"]),
            )
            super().__init__(
                generative_model=self.generative_model,
                amortizer=self.amortizer,
                configurator=generative_model.configurator,
                checkpoint_path=self.checkpoint_path + "/" + self.name,  # checkpoints,
                max_to_keep=self.trainer_setup["max_to_keep"],
                checkpoint_interval=self.trainer_setup["checkpoint_interval"],
                **kwargs,
            )
        else:
            super().__init__(generative_model=generative_model, **kwargs)

    def train_online(self, training_settings=None, save_file=False, **kwargs):
        if training_settings == "full":
            training_settings = FULL_TRAINING_SETTINGS
        if training_settings is None:
            eps = MyTrainer.train_settings["num_epochs"]
            its = MyTrainer.train_settings["iter_per_epoch"]
            bas = MyTrainer.train_settings["batch_size"]

            if "num_epochs" in kwargs:
                eps = training_settings["num_epochs"]
            if "iter_per_epoch" in kwargs:
                its = training_settings["iter_per_epoch"]
            if "batch_size" in kwargs:
                bas = training_settings["batch_size"]

        elif isinstance(training_settings, dict):
            eps = training_settings["num_epochs"]
            its = training_settings["iter_per_epoch"]
            bas = training_settings["batch_size"]

        loss = super().train_online(
            epochs=eps, iterations_per_epoch=its, batch_size=bas, **kwargs
        )
        if save_file:
            loss_plot = plot_losses(loss)
            plt.savefig(f"{self.figures_path}/loss_{self.name}.png")

        return loss
