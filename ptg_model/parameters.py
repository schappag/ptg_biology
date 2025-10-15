import numpy as np
from .utils import stim, sens, smooth_pw
from .model import releaseRate, rateAdj, DefaultParameters

def DefaultCalciumParams():
    return dict(lower=4.8, upper=5.2, unit="mg/dL")

def DefaultPhosphateParams():
    return dict(lower=3.2, upper=4.0, unit="mg/dL")

def DefaultVitaminDParams():
    return dict(lower=18, upper=61, unit="ng/L")

def SteadyState(C, Copt, Dopt):
    a = rateAdj(1, DefaultParameters("prolif"))
    ka = 0.001 * 60
    k2 = 0.03
    k1 = 4 * k2
    S2 = 1 / (1 + k2 / k1) * np.exp(-ka / a)
    S1 = k2 * S2 / k1
    S3 = rateAdj(1, DefaultParameters("prod")) * S1 / (
        releaseRate(C / 4, 1) + rateAdj(1, DefaultParameters("degrad"))
    )
    S4 = S3 * releaseRate(C / 4, 1) / rateAdj(1, DefaultParameters("clear"))
    return [S1, S2, S3, S4, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, Copt, Dopt, 0, 0, 0, 1]

def SteadyState_pat(C_pat, P_pat, D_pat, Copt, Popt, Dopt, PTH_pat, endpoints_D, endpoints_P):
    PTH = PTH_pat/9.434*3
    ka = 0.001*60
    k2SQ = 0.03
    k1SQ = 4*k2SQ
   
    kca = 0.5
    rca = 0.5
    kD = 0.001
    rD = 0.001
    D_pat= D_pat*smooth_pw(0, endpoints_D)
    P_pat = P_pat*smooth_pw(0, endpoints_P)
    aaCa = stim(C_pat-Copt, 'C')
    aaP  = stim(P_pat-Popt, 'P')
    aaD  = stim(D_pat-Dopt, 'D')

    yC = aaCa/(1+aaCa*np.sign(aaCa))
    yP = aaP/(1+aaP*np.sign(aaP))
    yD = aaD/(1+aaD*np.sign(aaD))

    a  = kca*(yC-2*yP)
    b  = kca*0.1
    aq = kD*(yD-2*yP)
    bq = kD*0.1

    ySC = (b-b*(bq-rD)/(aq-rD)-rca)/(a-b*bq/(aq-rD)-rca)
    ySD = (-rD+bq*(1-ySC))/(aq-rD)

    Csensed = sens(ySC, ySD)*C_pat
    Dsensed = sens(ySC, ySD)*D_pat

    aCa = stim(Csensed-Copt, 'C')

    Cstar  = aCa/(1+aCa*np.sign(aCa))
    Pstar  = aaP/(1+aaP*np.sign(aaP))
    CSstar = rca/(rca-50*kca*(Cstar-Pstar))

    ParamPNew = [Popt*0.323,0.3, 0.15, 4.5]
    KP = ParamPNew[0]
    aPhos = ParamPNew[1]
    bPhos = ParamPNew[2]
    gPhos = ParamPNew[3]

    FP = aPhos + (bPhos - aPhos) * (P_pat * 0.323) ** gPhos / (
        (P_pat * 0.323) ** gPhos + KP**gPhos
    )
    FP0 = aPhos + (bPhos - aPhos) * (Popt * 0.323) ** gPhos / (
        (Popt * 0.323) ** gPhos + KP**gPhos
    )

    rP = FP0 / FP
    GFR = np.exp(-10**(-3)*365*2)    
    S3 = PTH*rateAdj(GFR, DefaultParameters("clear"))/releaseRate(Csensed/4, rP)


    S1 = (releaseRate(Csensed/4, rP)+rateAdj(CSstar, DefaultParameters("degrad")))*S3 / rateAdj(CSstar, DefaultParameters("prod"))
    S2 = k1SQ*S1/k2SQ

    X = np.exp(ka/rateAdj(CSstar, DefaultParameters("prolif")))*(S1+S2)

    y_pat = np.array([S1, S2, S3, PTH, ySC, ySD, yC, yD, yP, CSstar, Cstar, CSstar, Cstar, CSstar, Cstar,
                      Csensed, Dsensed, Pstar, Pstar, Pstar, X])

    return y_pat