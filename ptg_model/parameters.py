'''
parameters.py
Utility functions for computing patient-specific steady states of the PTG model.
'''

import numpy as np
from ptg_model.core_functions import release_rate, rate_adj, defaultparameters
from ptg_model.utils import stim, sens, smooth_pw


def steady_state(c, copt, dopt):
    '''
    Compute the patient specific steady state of the PTG model given patient-specific calcium, phosphate, calcitriol, and PTH levels.
    Assumes phosphate and calcitriol are optimal.
    '''

    a = rate_adj(1, defaultparameters("prolif"))
    ka = 0.001 * 60
    k2sq = 0.03*60
    k1sq = 4 * k2sq
    s2 = 1 / (1 + k2sq / k1sq) * np.exp(-ka / a)
    s1 = k2sq * s2 / k1sq
    s3 = rate_adj(1, defaultparameters("prod")) * s1 / (
        release_rate(c / 4, 1) + rate_adj(1, defaultparameters("degrad"))
    )
    s4 = s3 * release_rate(c / 4, 1) / rate_adj(1, defaultparameters("clear"))
    return [s1, s2, s3, s4, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, copt, dopt, 0, 0, 0, 1]

def steadystate_pat(c_pat, p_pat, d_pat, copt, popt, dopt, pth_pat, endpoints_d, endpoints_p, gfr):
    '''
    Compute the patient specific steady state of the PTG model given patient-specific calcium, phosphate, calcitriol, and PTH levels.'''

    pth = pth_pat/9.434*3
    ka = 0.001*60
    k2sq = 0.03*60
    k1sq = 4*k2sq
   
    kca = 0.5
    rca = 0.5
    kd = 0.001
    rd = 0.001
    d_pat= d_pat*smooth_pw(0, endpoints_d)
    p_pat = p_pat*smooth_pw(0, endpoints_p)
    aaca = stim(c_pat-copt, 'c')
    aap  = stim(p_pat-popt, 'p')
    aad  = stim(d_pat-dopt, 'd')

    yc = aaca/(1+aaca*np.sign(aaca))
    yp = aap/(1+aap*np.sign(aap))
    yd = aad/(1+aad*np.sign(aad))

    a  = kca*(yc-2*yp)
    b  = kca*0.1
    aq = kd*(yd-2*yp)
    bq = kd*0.1

    ysc = (b-b*(bq-rd)/(aq-rd)-rca)/(a-b*bq/(aq-rd)-rca)
    ysd = (-rd+bq*(1-ysc))/(aq-rd)

    csensed = sens(ysc, ysd)*c_pat
    dsensed = sens(ysc, ysd)*d_pat

    aca = stim(csensed-copt, 'c')

    cstar  = aca/(1+aca*np.sign(aca))
    pstar  = aap/(1+aap*np.sign(aap))
    csstar = rca/(rca-50*kca*(cstar-pstar))

    param_phos = [popt*0.323,0.3, 0.15, 4.5]
    kp = param_phos[0]
    aphos = param_phos[1]
    bphos = param_phos[2]
    gphos = param_phos[3]

    fp = aphos + (bphos - aphos) * (p_pat * 0.323) ** gphos / (
        (p_pat * 0.323) ** gphos + kp**gphos
    )
    fp0 = aphos + (bphos - aphos) * (popt * 0.323) ** gphos / (
        (popt * 0.323) ** gphos + kp**gphos
    )

    rp = fp0 / fp
    s3 = pth*rate_adj(gfr, defaultparameters("clear"))/release_rate(csensed/4, rp)


    s1 = (release_rate(csensed/4, rp)+rate_adj(csstar, defaultparameters("degrad")))*s3 / rate_adj(csstar, defaultparameters("prod"))
    s2 = k1sq*s1/k2sq

    X = np.exp(ka/rate_adj(csstar, defaultparameters("prolif")))*(s1+s2)

    y_pat = np.array([s1, s2, s3, pth, ysc, ysd, yc, yd, yp, csstar, cstar, csstar, cstar, csstar, cstar,
                      csensed, dsensed, pstar, pstar, pstar, X])

    return y_pat