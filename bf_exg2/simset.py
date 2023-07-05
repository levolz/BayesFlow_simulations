import os
import numpy as np
from rpy2 import robjects


CHECKPOINT_DIR = "checkpoints"
FIGURES_DIR = "figures"
VAL_DIR = "validation_data"

if not os.path.exists(CHECKPOINT_DIR):
    os.makedirs(CHECKPOINT_DIR)

if not os.path.exists(FIGURES_DIR):
    os.makedirs(FIGURES_DIR)

if not os.path.exists(VAL_DIR):
    os.makedirs(VAL_DIR)

robjects.r(
    """
load("corey_data/corey_simulator.RData")
sds = diag(Sigma)
# load("corey_data/corey_indiv-fits_extracted.RData")
# means = corey_fits$log$means
# sds = corey_fits$log$sds
# sigma = corey_fits$log$cov
"""
)

global PRIOR_MEANS
global PRIOR_SDS
global PRIOR_SIGMA
global PRIOR_TRANS

global PRIOR_UPPER
global PRIOR_LOWER

global PARAM_NAMES

PRIOR_MEANS = np.asarray(robjects.r["mu"], order="F", dtype=np.float32)[:6]
PRIOR_SIGMA = np.asarray(robjects.r["Sigma"], order="F", dtype=np.float32)[:6, :6]
PRIOR_SDS = np.asarray(robjects.r["sds"], order="F", dtype=np.float32)[:6]

PRIOR_TRANS = np.asarray(robjects.r["te_final"], order="F", dtype=np.float32)[:6]

PRIOR_UPPER = np.asarray(robjects.r["tmv_upper"], order="F", dtype=np.float32)[:6]
PRIOR_LOWER = np.asarray(robjects.r["tmv_lower"], order="F", dtype=np.float32)[:6]

PARAM_NAMES = [
    r"$\mu_t$",
    r"$\sigma_t$",
    r"$\tau_t$",
    r"$\mu_s$",
    r"$\sigma_s$",
    r"$\tau_s$",
]
