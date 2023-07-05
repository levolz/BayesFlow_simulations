import numpy as np
from numba import njit


@njit(cache=True)
def exg(mu, sigma, tau):
    return np.random.normal(mu, sigma) + np.random.exponential(tau)


@njit
def rdiff(params):
    pass


@njit
def rt(params, SSD):
    ss_rt = SSD + exg(params[3], params[4], params[5])
    go_rt = exg(params[0], params[1], params[2])
    if go_rt < ss_rt:
        return go_rt
    else:
        return 0.0


@njit(cache=True)
def rts_stop(params, n_ss, SSD_start, SSD_step, SSD_ceil):
    ssd = SSD_start
    rts = np.zeros(n_ss, dtype=np.float32)
    SSD = np.zeros(n_ss, dtype=np.float32)
    SSD[0] = ssd
    for i in range(n_ss):
        SSD[i] = ssd
        rts[i] = rt(params, SSD[i])
        if rts[i] == 0.0:
            if ssd < SSD_ceil:
                ssd += SSD_step
        elif ssd > SSD_step:
            ssd -= SSD_step
    return rts, SSD


@njit(cache=True)
def rts(params, context):
    rts_ss, SSDs = rts_stop(
        params,
        context[1],
        context[3],
        context[4],
        context[5],
    )
    rts_go = np.asarray(
        [exg(params[0], params[1], params[2]) for _ in range(context[2])],
        dtype=np.float32,
    )
    RTs = np.concatenate((rts_ss, rts_go))
    SSDs = np.concatenate((SSDs, np.repeat(0.0, context[2])))
    return np.stack((RTs, SSDs), axis=1)


@njit(cache=True)
def simulate_batch(params, context):
    batch_size = params.shape[0]
    out = np.zeros((batch_size, context[0], 2), dtype=np.float32)
    for b in range(batch_size):
        out[b, :, :] = rts(params[b], context)
    return out


@njit(cache=True)
def uninf(params, context):
    batch_size = params.shape[0]
    s = np.zeros((batch_size, context[0]), dtype=np.float32)
    t = np.random.normal(0.0, 1.0, (batch_size, context[0]))
    return np.stack((t, s), axis=2)
