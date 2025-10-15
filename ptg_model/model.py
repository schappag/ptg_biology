import numpy as np
from .utils import stim, sens, smooth_pw

def rateAdj(C, parameterSet):
    r, A = parameterSet
    return np.where(C < 1, (r - A * r) * C + A * r, r)

def DefaultParameters(PTGFunction):
    if PTGFunction == "degrad":
        return [0.012 * 60, 0.1]
    if PTGFunction == "prolif":
        return [0.03 * 60, 1.2]
    if PTGFunction == "prod":
        return [6.6 / 0.1 * 60, 2]
    if PTGFunction == "clear":
        return [0.632 * 60, 0.2]

def releaseRate(C, rP):
    S_base, m, A, B = [1.21, 80, 0.14 * 60, 0.001 * 60]
    C_0 = 1.25
    S = S_base / 1.25 * C_0
    A *= rP
    return (A - B) / (1 + (C / S) ** m) + B

def deriv(t, y, endpoints_P, endpoints_D, Copt, Dopt, Popt,
          C_pat, P_pat, D_pat, S0, tf, GFR_in, y_pat):
    tauCa = 1
    tauD = 0.1
    tauP = 0.1
    kca = 0.5
    rca = 0.5
    kD = 0.001
    rD = 0.001
    ka = 0.001 * 60  # apoptosis rate

    k2 = 0.03
    k1 = 4 * k2
    tm = 24*30*12*3
    D = D_pat*smooth_pw(t/(tm), endpoints_D)
    C = C_pat*y[21]*y[22]
    P = P_pat*smooth_pw(t/(tm), endpoints_P)
    ParamPNew = [Popt*0.323,0.3, 0.15, 4.5]
    KP = ParamPNew[0]
    aPhos = ParamPNew[1]
    bPhos = ParamPNew[2]
    gPhos = ParamPNew[3]
    # Convert phosphate to mM
    FP = aPhos + (bPhos - aPhos) * (P * 0.323) ** gPhos / (
        (P * 0.323) ** gPhos + KP**gPhos
    )
    FP0 = aPhos + (bPhos - aPhos) * (Popt * 0.323) ** gPhos / (
        (Popt * 0.323) ** gPhos + KP**gPhos
    )

    rP = FP0 / FP

    dydt = np.zeros(23)

    dydt[0] = -k1 * y[0] + k2 * y[1]
    dydt[1] = (
        k1 * y[0]
        - k2 * y[1]
        - ka * y[1]
        + rateAdj(y[13], DefaultParameters("prolif"))
        * y[1]
        * np.log(y[20] / (y[1] + y[0]))
    )
    # releaseRate uses mmol/L for ionized calcium
    dydt[2] = (
        y[0] * rateAdj(y[11], DefaultParameters("prod"))
        - releaseRate(y[15] / 4, rP) * y[2]
        - rateAdj(y[9], DefaultParameters("degrad")) * y[2]
    )
    dydt[3] = releaseRate(y[15] / 4, rP) * y[2] - y[3] * rateAdj(
        GFR_in, DefaultParameters("clear")
    )
 

    dydt[4] = kca *((y[6] - 2*y[8])*y[4] + 0.1* (-1 + y[5])) + rca * (1 - y[4])

    dydt[5] = kD * ((y[7] - 2*y[8])*y[5] + 0.1* (-1 + y[4])) + rD * (1 - y[5])
    dydt[6] = (
        (stim(C- Copt, "C") * (1 - np.sign(stim(C- Copt, "C")) * y[6]) - y[6])
        * tauCa
        * 0.15
    )
    dydt[7] = (
        (stim(D - Dopt, "D") * (1 - np.sign(stim(D - Dopt, "D")) * y[7]) - y[7])
        * tauD
        * 0.15
    )
    dydt[8] = (
        (stim(P - Popt, "P") * (1 - np.sign(stim(P - Popt, "P")) * y[8]) - y[8])
        * tauP
        * 0.5
    )

    dydt[9] = 50 * kca * (y[10] - y[17]) * y[9] + rca * (1 - y[9])
    dydt[10] = (
        (
            stim(y[15] - Copt, "C") * (1 - np.sign(stim(y[15] - Copt, "C")) * y[10])
            - y[10]
        )
        * tauCa
        * 10
    )

    dydt[11] = 50 * kca * (y[12] - y[18]) * (y[11]) + rca * (1 - y[11])
    dydt[12] = (
        (
            stim(y[15] - Copt, "C") * (1 - np.sign(stim(y[15] - Copt, "C")) * y[12])
            - y[12]
        )
        * tauCa
        * 0.1
    )

    dydt[13] = 50 * kca * (y[14] - y[19]) * (y[13]) + rca * (1 - y[13])
    dydt[14] = (
        (
            stim(y[15] - Copt, "C") * (1 - np.sign(stim(y[15] - Copt, "C")) * y[14])
            - y[14]
        )
        * tauCa
        * 3.5
        * 10 ** (-1)
    )

    dydt[15] = sens(y[4], y[5]) * C - y[15]
    dydt[16] = sens(y[4], y[5]) * D - y[16]

    dydt[17] = (
        (stim(P - Popt, "P") * (1 - np.sign(stim(P - Popt, "P")) * y[17]) - y[17])
        * tauP
        * 1
    )
    dydt[18] = (
        (stim(P - Popt, "P") * (1 - np.sign(stim(P - Popt, "P")) * y[18]) - y[18])
        * tauP
        * 10 ** (-1)
        * 3
    )
    dydt[19] = (
        (stim(P - Popt, "P") * (1 - np.sign(stim(P - Popt, "P")) * y[19]) - y[19])
        * tauP
        * 10 ** (-1)*0.5
    )

    dydt[20] =  10 ** (-7) * (max(0, ((y[0] + y[1]) / S0 - 1)) ** (2 / 3))

    C_F = Copt / 4
    # pD = [0.001/4*C_F*0.4 ,0.001]
    # pDD = [0.05, 0.001]

    pD = [0.001/4*C_F*0.4 ,0.1]
    pDD = [0.05, 0.1]

    dydt[21] = pD[0] * (1 + pD[1] * np.tanh((y[3] - y_pat[3]) / y0[3]) * y[21])
    dydt[22] = pDD[0] * (1 + pDD[1] * np.tanh((D-D_pat - y_pat[3]) / y0[3]) * y[22])

    # dydt[21] = pD[0] * (1 + pD[1] * (y[3] - y_pat[3]) * y[21])
    #dydt[22] = pDD[0] * (1 + pDD[1] * (D-D_pat) -  y[22])
    threshold = 1e-6

    dydt[np.abs(dydt) < threshold] = 0
    return dydt
