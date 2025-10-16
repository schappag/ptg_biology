"""
model.py
System of ODEs describing PTG biology.
"""

import numpy as np
from ptg_model.utils import stim, sens, smooth_pw
from ptg_model.core_functions import rate_adj, defaultparameters, release_rate


def deriv(
    t,
    y,
    endpoints_p,
    endpoints_d,
    copt,
    dopt,
    popt,
    c_pat,
    p_pat,
    d_pat,
    s0,
    tm,
    gfr_in,
    y_pat,
    calcium_clamp=True,
):
    """
    Defines the system of ODEs describing PTG biology.

    Parameters
    ----------
    *args : tuple
        Model parameters including phosphate, calcitriol, and calcium input.

    Returns
    -------
    function
        A system of ODEs suitable for `solve_ivp`.
    """
    tauca = 1  # Time parameters
    taud = 0.1
    taup = 0.1
    kca = 0.5
    rca = 0.5
    kd = 0.001
    rd = 0.001

    ka = 0.001 * 60  # apoptosis rate

    k2 = 0.03 * 60  # Transition rates
    k1 = 4 * k2
    d = d_pat * smooth_pw(t / (tm), endpoints_d)
    if calcium_clamp:
        c = c_pat
    else:
        c = c_pat * y[21] * y[22]
    p = p_pat * smooth_pw(t / (tm), endpoints_p)

    parameter_phos = [popt * 0.323, 0.3, 0.15, 4.5]
    kp = parameter_phos[0]
    aphos = parameter_phos[1]
    bphos = parameter_phos[2]
    gphos = parameter_phos[3]

    # convert phosphate to mM
    fp = aphos + (bphos - aphos) * (p * 0.323) ** gphos / (
        (p * 0.323) ** gphos + kp**gphos
    )
    fp0 = aphos + (bphos - aphos) * (popt * 0.323) ** gphos / (
        (popt * 0.323) ** gphos + kp**gphos
    )

    rp = fp0 / fp

    dydt = np.zeros(23)

    dydt[0] = -k1 * y[0] + k2 * y[1]
    dydt[1] = (
        k1 * y[0]
        - k2 * y[1]
        - ka * y[1]
        + rate_adj(y[13], defaultparameters("prolif"))
        * y[1]
        * np.log(y[20] / (y[1] + y[0]))
    )
    # release_rate uses mmol/L for ionized calcium
    dydt[2] = (
        y[0] * rate_adj(y[11], defaultparameters("prod"))
        - release_rate(y[15] / 4, rp) * y[2]
        - rate_adj(y[9], defaultparameters("degrad")) * y[2]
    )
    dydt[3] = release_rate(y[15] / 4, rp) * y[2] - y[3] * rate_adj(
        gfr_in, defaultparameters("clear")
    )

    dydt[4] = kca * ((y[6] - 2 * y[8]) * y[4] + 0.1 * (-1 + y[5])) + rca * (1 - y[4])

    dydt[5] = kd * ((y[7] - 2 * y[8]) * y[5] + 0.1 * (-1 + y[4])) + rd * (1 - y[5])
    dydt[6] = (
        (stim(c - copt, "c") * (1 - np.sign(stim(c - copt, "c")) * y[6]) - y[6])
        * tauca
        * 0.15
    )
    dydt[7] = (
        (stim(d - dopt, "d") * (1 - np.sign(stim(d - dopt, "d")) * y[7]) - y[7])
        * taud
        * 0.15
    )
    dydt[8] = (
        (stim(p - popt, "p") * (1 - np.sign(stim(p - popt, "p")) * y[8]) - y[8])
        * taup
        * 0.5
    )

    dydt[9] = 50 * kca * (y[10] - y[17]) * y[9] + rca * (1 - y[9])
    dydt[10] = (
        (
            stim(y[15] - copt, "c") * (1 - np.sign(stim(y[15] - copt, "c")) * y[10])
            - y[10]
        )
        * tauca
        * 10
    )

    dydt[11] = 50 * kca * (y[12] - y[18]) * (y[11]) + rca * (1 - y[11])
    dydt[12] = (
        (
            stim(y[15] - copt, "c") * (1 - np.sign(stim(y[15] - copt, "c")) * y[12])
            - y[12]
        )
        * tauca
        * 0.1
    )

    dydt[13] = 50 * kca * (y[14] - y[19]) * (y[13]) + rca * (1 - y[13])
    dydt[14] = (
        (
            stim(y[15] - copt, "c") * (1 - np.sign(stim(y[15] - copt, "c")) * y[14])
            - y[14]
        )
        * tauca
        * 3.5
        * 10 ** (-1)
    )

    dydt[15] = sens(y[4], y[5]) * c - y[15]
    dydt[16] = sens(y[4], y[5]) * d - y[16]

    dydt[17] = (
        (stim(p - popt, "p") * (1 - np.sign(stim(p - popt, "p")) * y[17]) - y[17])
        * taup
        * 1
    )
    dydt[18] = (
        (stim(p - popt, "p") * (1 - np.sign(stim(p - popt, "p")) * y[18]) - y[18])
        * taup
        * 10 ** (-1)
        * 3
    )
    dydt[19] = (
        (stim(p - popt, "p") * (1 - np.sign(stim(p - popt, "p")) * y[19]) - y[19])
        * taup
        * 10 ** (-1)
        * 0.5
    )

    dydt[20] = 10 ** (-5) * (max(0, ((y[0] + y[1]) / s0 - 1)) ** (2 / 3))

    c_f = copt / 4

    pd = [0.001 / 4 * c_f * 0.4, 0.1]
    pdd = [0.05, 0.1]

    target21 = 1 + np.tanh(pd[1] * (y[3] - y_pat[3]))
    target22 = 1 + np.tanh(pdd[1] * (d - d_pat))

    dydt[21] = pd[0] * (target21 - y[21])
    dydt[22] = pdd[0] * (target22 - y[22])

    threshold = 1e-12

    dydt[np.abs(dydt) < threshold] = 0
    return dydt
