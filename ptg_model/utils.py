import numpy as np

def smooth_pw(x, endpoints, alpha=80):
    endpoints_x, endpoints_y = endpoints
    beta = endpoints_x[1:-1]
    JP = (endpoints_y[1:] - endpoints_y[:-1]) / (endpoints_x[1:] - endpoints_x[:-1])
    bP = (JP[-1] + JP[0]) / 2
    cP = (JP[1:] - JP[:-1]) / 2
    aP = endpoints_y[0] - np.sum(cP * np.abs(beta))
    AP = aP - np.sum(cP * beta)
    BP = bP + np.sum(cP)
    return AP + BP * x + np.sum(cP * np.log(1 + np.exp(-alpha * (x - beta)))) * 2 / alpha

def smooth_pw_matrix(x, endpoints, alpha=100):
    return np.array([smooth_pw(val, endpoints, alpha) for val in x])

def stim(Val, Param):
    if Param == "C":
        C1, C2, k, L = -2.2, 2.2, 3, 1
    elif Param == "P":
        C1, C2, k, L = -2.5, 2.5, 2.5, 1
    elif Param == "D":
        C1, C2, k, L = -30, 30, 0.1, 1
    s = L/(1+np.exp(-k*(Val-C1))) + L/(1+np.exp(-k*(Val-C2))) - L
    cutoff = L/(1+np.exp(-k*(C2-C1)/2)) + L/(1+np.exp(-k*(C1-C2)/2)) - L
    return np.where(np.abs(s) > cutoff, s, 0)

def sens(C, D):
    endpoints = np.array([[0, 0.5, 1, 2, 10],
                          [0.65, 0.7, 1, 1.01, 1.05]])
    avg = (C + D) / 2
    return smooth_pw(avg, endpoints) / smooth_pw(1, endpoints)
