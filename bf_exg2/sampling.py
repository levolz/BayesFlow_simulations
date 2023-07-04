import numpy as np
from numba import njit
from numba_stats import norm


@njit(parallel=True)
def mvn(mu, sigma, size):
    """
    Multivariate normal distribution.

    Parameters
    ----------
    mu : array_like
        Mean of the distribution.
    sigma : array_like
        Covariance matrix of the distribution.

    Returns
    -------
    out : ndarray
        Random samples from the multivariate normal distribution.
    """
    n = mu.shape[0]

    z = np.zeros((size, n), dtype=np.float32)
    for i in range(n):
        z[:, i] = np.random.normal(0.0, 1.0, size=size)
    L = np.asarray(np.linalg.cholesky(sigma), dtype=np.float32)
    return mu + np.dot(z, L.T)


@njit(parallel=True)
def uvn(mu, sd, size):
    n = mu.shape[0]
    out = np.zeros((size, n), dtype=np.float32)
    for i in range(n):
        out[:, i] = np.random.normal(mu[i], sd[i], size=size)
    return out


@njit(parallel=True)
def lhs(k, size):
    intervals = np.linspace(0.0, 1.0, size + 1)[:-1]
    spacing = intervals[1]
    out = np.zeros((size, k), dtype=np.float32)
    for j in range(k):
        out[:, j] = intervals + spacing * np.random.rand(size)
        np.random.shuffle(out[:, j])
    return out


@njit(parallel=True)
def uvu(lower, upper, size):
    n = lower.shape[0]
    out = np.zeros((size, n), dtype=np.float32)
    for i in range(n):
        l = lower[i]
        u = upper[i]
        if l <= u:
            out[:, i] = np.random.uniform(l, u, size=size)
        else:
            out[:, i] = np.random.uniform(u, l, size=size)
    return out


@njit(parallel=True)
def uvu_lhs(lower, upper, size):
    k = lower.shape[0]
    intervals = np.linspace(0.0, 1.0, size + 1)[:-1]
    spacing = intervals[1]

    out = np.zeros((size, k), dtype=np.float32)
    for j in range(k):
        l = lower[j]
        u = upper[j]
        samp = intervals + spacing * np.random.rand(size)
        if l <= u:
            out[:, j] = l + samp * (u - l)
        else:
            out[:, j] = u + samp * (l - u)
        np.random.shuffle(out[:, j])
    return out


@njit(parallel=True)
def uvn_lhs(mu, sd, size):
    k = mu.shape[0]
    intervals = np.linspace(0.0, 1.0, size + 1)[:-1]
    spacing = intervals[1]

    out = np.zeros((size, k), dtype=np.float32)
    for j in range(k):
        out[:, j] = norm.ppf(intervals + spacing * np.random.rand(size), 0.0, 1.0)
        np.random.shuffle(out[:, j])
    return mu + sd * out


@njit(parallel=True)
def mvn_lhs(mu, sigma, size):
    k = mu.shape[0]
    intervals = np.linspace(0.0, 1.0, size + 1)[:-1]
    spacing = intervals[1]

    out = np.zeros((size, k), dtype=np.float32)
    for j in range(k):
        out[:, j] = norm.ppf(intervals + spacing * np.random.rand(size), 0.0, 1.0)
        np.random.shuffle(out[:, j])

    L = np.asarray(np.linalg.cholesky(sigma), dtype=np.float32)
    return mu + np.dot(out, L.T)
