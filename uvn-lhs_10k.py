from bf_exg2.exg2 import MyPrior, MySimulator, MyGenerativeModel, MyTrainer
from bf_exg2.diagnostics import Diagnostics

import pickle as pkl

SIZE = 10000
NAME = "uvn-lhs_10k"


# Path: bf_exg2\tests\generate.py
def main():
    prior = MyPrior("uvn", lhs_sampling=True)
    simulator = MySimulator(n_trials=SIZE)
    model = MyGenerativeModel(prior, simulator)
    trainer = MyTrainer(model, name=NAME)
    loss = trainer.train_online(training_settings=None, save_file=True)
    pkl.dump(loss, open(f"{NAME}.pkl", "wb"))
    diagnostics = Diagnostics(trainer)
    diagnostics.validation_sample(100, 100)
    diagnostics.plot_posterior_recovery(save=True)


if __name__ == "__main__":
    main()
