from bf_exg2.exg2 import MyPrior, MySimulator, MyGenerativeModel, MyTrainer
from bf_exg2.diagnostics import Diagnostics, Inference

import pickle as pkl

TEST_TRAINING_SETTINGS = {
    "num_epochs": 2,
    "iter_per_epoch": 8,
    "batch_size": 32,
}

SIZE = {
    "low": 200,
    "med": 600,
    "high": 1000,
}

NAME = "test"


# Path: bf_exg2\tests\generate.py
def main():
    prior = MyPrior("uvu", lhs_sampling=False, trans=True)
    simulator = MySimulator(n_trials=(80, 1200), uninformative=False)
    model = MyGenerativeModel(prior, simulator)
    trainer = MyTrainer(model, name=NAME)
    loss = trainer.train_online(
        training_settings=TEST_TRAINING_SETTINGS, save_file=True
    )
    pkl.dump(loss, open(f"{NAME}.pkl", "wb"))
    diagnostics = Diagnostics(trainer)
    diagnostics.validation_sample(2, 2, SIZE)
    diagnostics.plot_posterior_recovery(save=True)


if __name__ == "__main__":
    main()
