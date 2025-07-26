# Simulation project for stop-signal paradigm in BayesFlow
See [`poster.pdf`](https://github.com/levolz/BayesFlow_simulations/blob/main/poster.pdf) for a summary of results.

`bf_exg2` contains the packaged simulation and analysis functions, extending on [BayesFlow](https://bayesflow.org/) functionality in `bayesflow` (2022 version)
`DMC` is a copy of the functionality introduced in [Heathcote et al. (2019)](https://doi.org/10.3758/s13428-018-1067-y) --- `.cpp` files contain self-written compiled code to accelerate computations
Files in `main` each run a specific exGaussian model with given simulation settings
