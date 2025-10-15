from scipy.integrate import solve_ivp

def run_simulation(deriv, y0, t_span, *, args=(), method="BDF"):
    return solve_ivp(deriv, t_span, y0, args=args, method=method)
