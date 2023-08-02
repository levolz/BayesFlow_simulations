from bf_exg2.exg2 import MyPrior, MySimulator, MyGenerativeModel, MyTrainer
from bf_exg2.diagnostics import Diagnostics, Inference

import pickle as pkl
import numpy as np
import os
import csv

SIZE = {
    "low": 200,
    "med": 600,
    "high": 1000,
}

NAME = "uvn-lhs-12"


# Path: bf_exg2\tests\generate.py
def main():
    prior = MyPrior("uvn", lhs_sampling=True)
    simulator = MySimulator(n_trials=(80, 1200), uninformative=False)
    model = MyGenerativeModel(prior, simulator)
    trainer = MyTrainer(model, name=NAME)
    # loss = trainer.train_online(training_settings="full", save_file=True)
    # pkl.dump(loss, open(f"{NAME}.pkl", "wb"))
    diagnostics = Diagnostics(trainer)
    diagnostics.validation_sample(100, 100, SIZE)
    diagnostics.plot_posterior_recovery(save=True)
    inference = Inference(trainer, name="comparison")
    samples = np.zeros((114, 100, 6))
    for i, file in enumerate(os.listdir("corey_data/bf_new")):
        with open(f"corey_data/bf_new/{file}", "r") as f:
            reader = csv.reader(f)
            data = np.asarray([list(reader)[1:]], dtype=np.float32)
        samples[i, :, :] = inference.sample_posterior(
            n_samples=100, data=data, parameters=None, _return=True
        )

    with open(f"corey_data/uvn-dmc_params.csv", "r") as f:
        reader = csv.reader(f)
        tmp = list()
        next(reader)
        for row in reader:
            tmp.append(row)
        dmc_params = np.asarray(tmp, dtype=np.float32)

    diagnostics.plot_posterior_recovery(
        post_samples=samples * model.prior_sd + model.prior_mu,
        validation_sims={"parameters": dmc_params},
        uncertainty_agg=None,
        save=True,
        standardized=True,
        x_label="DMC",
        y_label="BayesFlow",
    )


if __name__ == "__main__":
    main()
