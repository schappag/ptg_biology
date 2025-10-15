"""
core_functions.py
-----------------
Core PTG model helper functions shared between model.py and parameters.py.

Includes:
- rate_adj: adjusts PTG parameters based on input concentration.
- defaultparameters: provides default kinetic parameter sets.
- release_rate: defines the sigmoidal PTH release function.
"""

import numpy as np


def rate_adj(c, parameterset):
    """Adjust PTG parameters based on calcium/phosphate input.

    Parameters
    ----------
    c : float or ndarray
        Current input concentration (e.g., calcium).
    parameterset : list or tuple
        Two-element list containing the base rate and adjustment factor.

    Returns
    -------
    float or ndarray
        Adjusted parameter value.
    """
    r, a = parameterset
    return np.where(c < 1, (r - a * r) * c + a * r, r)


def defaultparameters(ptgfunction):
    """Return default parameters for the four PTG functions.

    Parameters
    ----------
    ptgfunction : str
        One of {'degrad', 'prolif', 'prod', 'clear'}.

    Returns
    -------
    list
        Default parameter pair [rate, adjustment].
    """
    if ptgfunction == "degrad":
        return [0.012 * 60, 0.1]
    if ptgfunction == "prolif":
        return [0.03 * 60, 2]
    if ptgfunction == "prod":
        return [6.6 / 0.1 * 60, 2]
    if ptgfunction == "clear":
        return [0.632 * 60, 0.2]
    raise ValueError(f"Unknown PTG function: {ptgfunction}")


def release_rate(c, rp, copt=1.25):
    """Sigmoidal PTH release function.

    Parameters
    ----------
    c : float
        Calcium concentration (mmol/L).
    rp : float
        Scaling factor from phosphate regulation.
    copt : float, optional
        Optimal calcium value. Default is 1.25.

    Returns
    -------
    float
        Release rate value.
    """
    s_base, m, a, b = [1.22, 100, 0.14 * 60, 0.001 * 60]
    s = s_base / 1.25 * copt
    a *= rp
    return (a - b) / (1 + (c / s) ** m) + b
