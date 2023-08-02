from bf_exg2.exg2 import MyPrior, MySimulator, MyGenerativeModel, MyTrainer
from bf_exg2.diagnostics import Diagnostics, Inference

import pickle as pkl

SIZE = 10000
NAME = "uvn-lhs_10k-12"


# Path: bf_exg2\tests\generate.py
def main():
    prior = MyPrior("uvn", lhs_sampling=True)
    simulator = MySimulator(n_trials=SIZE)
    model = MyGenerativeModel(prior, simulator)
    trainer = MyTrainer(model, name=NAME)
    # loss = trainer.train_online(training_settings=None, save_file=True)
    # pkl.dump(loss, open(f"{NAME}.pkl", "wb"))
    diagnostics = Diagnostics(trainer)
    diagnostics.validation_sample(100, 100)
    diagnostics.plot_posterior_recovery(save=True)

    inference = Inference(trainer)
    inference.sample(100)
    inference.save_for_dmc()
    inference.sample_posterior(100)
    pkl.dump(inference.posterior_samples, open(f"{NAME}_post-samples.pkl", "wb"))


if __name__ == "__main__":
    main()
