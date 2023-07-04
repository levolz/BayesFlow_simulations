import importlib.resources as pkg_res
from bf_exg2 import data
import pickle

PARAM_NAMES = [
    r"$\mu_t$",
    r"$\sigma_t$",
    r"$\tau_t$",
    r"$\mu_s$",
    r"$\sigma_s$",
    r"$\tau_s$",
]

DIRS = {
    "checkpoints": "checkpoints",
    "figures": "figures",
    "validation_data": "validation_data",
}

STAIRCASE = {
    "SSD_start": 0.25,
    "SSD_step": 0.05,
    "SSD_ceil": 1.0,
    "prop_ss": 0.25,
}

TRAINER_SETUP = {
    "summary_dim": 4,
    "num_coupling_layers": 4,
    "max_to_keep": 8,
    "checkpoint_interval": 32,
}

TRAINING_SETTINGS = {
    "num_epochs": 128,
    "iter_per_epoch": 256,
    "batch_size": 64,
}


FULL_TRAINING_SETTINGS = {
    "num_epochs": 1024,
    "iter_per_epoch": 1024,
    "batch_size": 64,
}


def default_prior():
    file = pkg_res.files(data) / "prior.pkl"
    with file.open("rb") as f:
        dat = pickle.load(f)
    return dat
