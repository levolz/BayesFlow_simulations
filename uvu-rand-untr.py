from bf_exg2.exg2 import MyPrior, MySimulator, MyGenerativeModel, MyTrainer
from bf_exg2.diagnostics import Diagnostics

import pickle as pkl

SIZE = {
    "low": 200,
    "med": 600,
    "high": 1000,
}

NAME = "uvu-rand-untr"


# Path: bf_exg2\tests\generate.py
def main():
    prior = MyPrior("uvu", lhs_sampling=True, trans=False)
    simulator = MySimulator(n_trials=(80, 1200), uninformative=False)
    model = MyGenerativeModel(prior, simulator)
    trainer = MyTrainer(model, name=NAME)
    loss = trainer.train_online(training_settings="full", save_file=True)
    pkl.dump(loss, open(f"{NAME}.pkl", "wb"))
    diagnostics = Diagnostics(trainer)
    diagnostics.validation_sample(100, 100, SIZE)
    diagnostics.plot_posterior_recovery(save=True)


if __name__ == "__main__":
    main()
