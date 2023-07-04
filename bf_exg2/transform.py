import numpy as np
from numba import njit


@njit(parallel=True)
def boxcox(x, trans):
    """
    Box-Cox transformation of x

    Parameters
    ----------
    x : np.ndarray (n, k)
        Data to be transformed
    trans : np.ndarray (k,)
        Transformation parameters
    """
    out = np.zeros_like(x)
    for j, t in enumerate(trans):
        if t == 0.0:
            out[:, j] = np.log(x[:, j])
        else:
            out[:, j] = x[:, j] ** t
    return out


@njit(parallel=True)
def boxcox_inv(x, trans):
    """
    Inverse Box-Cox transformation of x

    Parameters
    ----------
    x : np.ndarray (n, k)
        Data to be transformed
    trans : np.ndarray (k,)
        Transformation parameters
    """
    out = np.zeros_like(x)
    for j, t in enumerate(trans):
        if t == 0.0:
            out[:, j] = np.exp(x[:, j])
        else:
            out[:, j] = x[:, j] ** (1 / t)
    return out


@njit(parallel=True)
def windsorise(x, lower, upper):
    x = np.minimum(x, upper)
    x = np.maximum(x, lower)
    return x


@njit(parallel=True)
def retransform(x, bounds, trans):
    return boxcox_inv(windsorise(x, bounds[0], bounds[1]), trans)
