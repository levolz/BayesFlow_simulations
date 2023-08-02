import numpy as np
from rpy2 import robjects
import os

import pickle

print(os.getcwd())

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

PRIOR = {
    "mu": np.asarray(robjects.r["mu"], order="F", dtype=np.float32)[:6],
    "sd": np.sqrt(np.asarray(robjects.r["sds"], order="F", dtype=np.float32)[:6]),
    "sigma": np.asarray(robjects.r["Sigma"], order="F", dtype=np.float32)[:6, :6],
    "trans": np.asarray(robjects.r["te_final"], order="F", dtype=np.float32)[:6],
    "upper": np.asarray([1e6, 0, 1, 1, 1e-6, 0], order="F", dtype=np.float32),
    # np.asarray(robjects.r["tmv_upper"], order="F", dtype=np.float32)[:6],
    "lower": np.asarray([1, -1e6, 1e-6, 1e-6, -1e6, -1e6], order="F", dtype=np.float32),
    # np.asarray(robjects.r["tmv_lower"], order="F", dtype=np.float32)[:6],
}
with open(
    "bf_exg2/data/prior.pkl",
    "wb",
) as f:
    pickle.dump(PRIOR, f)
