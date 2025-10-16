"""
ptg_model.utils
---------------
Utility functions for smooth piecewise interpolation, stimulation, and sensitivity
calculations used in the pTG (parathyroid gland) model.
"""

import numpy as np


def smooth_pw(x, endpoints, alpha=80):
    """
    compute a smooth piecewise interpolation between defined endpoints.

    This function uses a smooth logistic approximation to create continuous,
    differentiable transitions between line segments defined by (x, y) endpoints.

    parameters
    ----------
    x : float or ndarray
        Input value(s) where the smooth function should be evaluated.
    endpoints : ndarray
        Array of shape (2, N) where the first row contains the x-coordinates
        and the second row contains the y-coordinates of the endpoints.
    alpha : float, optional
        Smoothness parameter controlling the sharpness of transitions.
        Higher alpha → steeper transitions. Default is 80.

    Returns
    -------
    float or ndarray
        Smoothed function value(s) corresponding to `x`.
    """
    endpoints_x, endpoints_y = endpoints
    beta = endpoints_x[1:-1]
    jp = (endpoints_y[1:] - endpoints_y[:-1]) / (endpoints_x[1:] - endpoints_x[:-1])
    bp = (jp[-1] + jp[0]) / 2
    cp = (jp[1:] - jp[:-1]) / 2
    ap = endpoints_y[0] - np.sum(cp * np.abs(beta))
    a_p = ap - np.sum(cp * beta)
    b_p = bp + np.sum(cp)
    z = -alpha * (x - beta)
    return a_p + b_p * x + np.sum(cp * np.log1p(np.exp(z))) * 2 / alpha


def smooth_pw_matrix(x, endpoints, alpha=100):
    """
    Vectorized version of `smooth_pw`.

    Evaluates the smooth piecewise interpolation for an array of input values.

    Parameters
    ----------
    x : array_like
        Input values at which to evaluate the smooth function.
    endpoints : ndarray
        Array of shape (2, N) containing x and y endpoints.
    alpha : float, optional
        Smoothness parameter (default = 100).

    Returns
    -------
    ndarray
        Smoothed function values at all points in `x`.
    """
    return np.array([smooth_pw(val, endpoints, alpha) for val in x])


def stim(val, param):
    """
    compute the stimulation function for calcium, phosphate, or calcitriol.

    parameters
    ----------
    val : float or ndarray
        Input variable value(s) (e.g., calcium concentration).
    param : {'c', 'p', 'Dc'}
        Type of parameter to simulate:
        - 'c' : calcium
        - 'p' : phosphate
        - 'c' : calcitriol

    Returns
    -------
    ndarray
        Stimulation values with cutoff applied to small responses.
    """
    if param == "c":
        c1, c2, k, l = -2.2, 2.2, 3, 1
    elif param == "p":
        c1, c2, k, l = -2.5, 2.5, 2.5, 1
    elif param == "d":
        c1, c2, k, l = -30, 30, 0.1, 1
    else:
        c1, c2, k, l = 0.0, 0.0, 1.0, 1.0  # Default fallback values

    s = l / (1 + np.exp(-k * (val - c1))) + l / (1 + np.exp(-k * (val - c2))) - l
    cutoff = (
        l / (1 + np.exp(-k * ((c2 - c1) / 2)))
        + l / (1 + np.exp(-k * ((c1 - c2) / 2)))
        - l
    )
    return np.where(np.abs(s) > cutoff, s, 0)


def sens(c, d):
    """
    compute a sensitivity scaling factor

    parameters
    ----------
    c : float or ndarray
        calcium-related value(s).
    D : float or ndarray
        calcitriol-related value(s).

    Returns
    -------
    float or ndarray
        Normalized sensitivity scaling factor.
    """
    endpoints = np.array([[0, 0.5, 1, 2, 10], [0.65, 0.7, 1, 1.01, 1.05]])
    avg = (c + d) / 2
    return smooth_pw(avg, endpoints) / smooth_pw(1, endpoints)
