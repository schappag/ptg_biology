import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer

from config import simulation_config as cfg
from ptg_model.utils import smooth_pw
from ptg_model.parameters import SteadyState, SteadyState_pat
from ptg_model.model import deriv
from ptg_model.simulation import run_simulation
from ptg_model.plot import plot_iPTH, plot_states


def get_params_interactive(cfg):
    """Ask user for overrides, keep defaults if Enter is pressed."""
    def ask(key, default):
        value = input(f"{key} (default {default}): ")
        if value == "":
            return default
        try:
            return type(default)(value)
        except ValueError:
            return value

    cfg["unit_Ca"] = ask("Unit of Calcium (mg/dL, mmol/L, mEq/L)", cfg["unit_Ca"])
    cfg["lower_Ca"] = ask("Lower Ca ref", cfg["lower_Ca"])
    cfg["upper_Ca"] = ask("Upper Ca ref", cfg["upper_Ca"])
    cfg["measurement_type"] = ask("Calcium type (total/ionized)", cfg["measurement_type"])

    cfg["unit_P"] = ask("Unit of Phosphate (mg/dL, mmol/L)", cfg["unit_P"])
    cfg["lower_P"] = ask("Lower P ref", cfg["lower_P"])
    cfg["upper_P"] = ask("Upper P ref", cfg["upper_P"])

    cfg["unit_D"] = ask("Unit of Vitamin D (ng/L, pmol/L)", cfg["unit_D"])
    cfg["lower_D"] = ask("Lower D ref", cfg["lower_D"])
    cfg["upper_D"] = ask("Upper D ref", cfg["upper_D"])

    cfg["unit_iPTH"] = ask("Unit of iPTH (pg/mL, pmol/L)", cfg["unit_iPTH"])
    cfg["vintage"] = ask("Dialysis vintage (years)", cfg["vintage"])
    cfg["P_control"] = ask("Phosphate control (poor/good/excellent)", cfg["P_control"])
    cfg["D_control"] = ask("Calcitriol therapy (yes/no)", cfg["D_control"])

    cfg["tf"] = ask("Simulation horizon", cfg["tf"])
    cfg["tf_unit"] = ask("Horizon unit (h/d/m)", cfg["tf_unit"])

    cfg["target_P"] = ask("Target phosphate", cfg["target_P"])
    cfg["target_D"] = ask("Target calcitriol", cfg["target_D"])
    cfg["time_to_target_P"] = ask("Time to target P (months)", cfg["time_to_target_P"])
    cfg["time_to_target_D"] = ask("Time to target D (months)", cfg["time_to_target_D"])

    return cfg


def main():
    """
    Main function to set up and run the PTG biology simulation.
    This function interacts with the user to optionally override default simulation parameters,
    performs unit conversions and parameter calculations for calcium (Ca), phosphate (P), and vitamin D (D),
    sets up simulation endpoints and initial conditions, and runs the simulation using the specified ODE solver.
    After the simulation, it generates plots for iPTH and state variables.
    Steps performed:
    1. Prompts user to override default parameters or use existing configuration.
    2. Extracts and processes simulation parameters, including unit conversions and control settings.
    3. Calculates optimal and patient-specific values for Ca, P, and D.
    4. Sets up time normalization and target endpoints for P and D.
    5. Initializes the system's steady state and patient state.
    6. Prepares parameters and runs the simulation using the `run_simulation` function.
    7. Plots the results for iPTH and state variables.
    Assumes the existence of the following functions and variables:
    - `cfg`, `get_params_interactive`, `SteadyState`, `SteadyState_pat`, `deriv`, `run_simulation`
    - `plot_iPTH`, `plot_states`, `plt`, `np`, `timer`
    """
    # Ask whether to override
    override = input("Do you want to override default parameters? (y/n): ").strip().lower()
    if override == "y":
        params = get_params_interactive(cfg.copy())
    else:
        params = cfg.copy()

    # === Simulation parameter setup ===
    unit_Ca = params["unit_Ca"]
    lower_Ca = params["lower_Ca"]
    upper_Ca = params["upper_Ca"]
    measurement_type = params["measurement_type"]

    mean_Ca = (lower_Ca + upper_Ca) / 2
    if measurement_type.lower() == "total":
        mean_Ca *= 0.5
    if unit_Ca == "mmol/L":
        Copt = mean_Ca * 4
    elif unit_Ca == "mEq/L":
        Copt = mean_Ca * 2
    else:
        Copt = mean_Ca

    unit_P = params["unit_P"]
    lower_P = params["lower_P"]
    upper_P = params["upper_P"]
    mean_P = (lower_P + upper_P) / 2
    if unit_P == "mmol/L":
        Popt = mean_P / 0.323
    else:
        Popt = mean_P

    unit_D = params["unit_D"]
    lower_D = params["lower_D"]
    upper_D = params["upper_D"]
    mean_D = (lower_D + upper_D) / 2
    if unit_D == "pmol/L":
        Dopt = mean_D / 2.4001
    else:
        Dopt = mean_D

    unit_iPTH = params["unit_iPTH"]
    if unit_iPTH == "pmol/L":
        iPTH_factor = 0.106
    else:
        iPTH_factor = 1

    vintage = params["vintage"]
    P_control = params["P_control"]
    D_control = params["D_control"]

    if P_control == "excellent":
        P_pat = 4.6
    elif P_control == "good":
        P_pat = 5.5
    else:
        P_pat = 8

    if D_control == "no":
        D_pat = 0.5 * lower_D
    else:
        D_pat = 0.8 * Dopt

    tf = params["tf"]
    tf_unit = params["tf_unit"]
    if tf_unit == "d":
        tf *= 24
    elif tf_unit == "m":
        tf *= 24 * 30

    target_P = params["target_P"]
    target_D = params["target_D"]
    time_to_target_P = params["time_to_target_P"] * 24 * 30
    time_to_target_D = params["time_to_target_D"] * 24 * 30

    if unit_P == "mmol/L":
        target_P /= 0.323
    if unit_D == "pmol/L":
        target_D /= 2.4001

    t_norm_P = float(24*30*12*2)
    if tf < time_to_target_P:
        t_norm_P = time_to_target_P + 1

    t_norm_D = float(24*30*12*2)
    if tf < time_to_target_D:
        t_norm_D = time_to_target_D + 1

    C_pat = Copt * 0.95

    endpoints_P = np.array(
        [[0, time_to_target_P / t_norm_P, 1],
         [1, target_P / P_pat, target_P / P_pat]]
    )
    endpoints_D = np.array(
        [[0, time_to_target_D / t_norm_D, 1],
         [1, target_D / D_pat, target_D / D_pat]]
    )

    y0 = SteadyState(Copt, Copt, Dopt)
    PTHopt = y0[3] / 3 * 9.434
    S0 = y0[0] + y0[1]

    PTH_pat = 900
    y_pat = SteadyState_pat(C_pat, P_pat, D_pat, Copt, Popt, Dopt, PTH_pat,
                            endpoints_D, endpoints_P)
    y_pat = np.append(y_pat, [1, 1])

    GFR_in = np.exp(-10**(-3)*365*2)
    p = (endpoints_P, endpoints_D, Copt, Dopt, Popt,
         C_pat, P_pat, D_pat, S0, tf, GFR_in, y_pat)

    d = deriv(0, y_pat, *p)
    print(f'deriv:{d}')

    # === Simulation ===
    t_span = (0, tf)
    start = timer()
    sn = run_simulation(deriv, y_pat, t_span, args=p, method="BDF")
    end = timer()

    # === Plots ===
    plot_iPTH(sn, unit_iPTH, abs_rel="abs")
    plot_states(sn, indices=range(0, 4))
    plt.show()


if __name__ == "__main__":
    main()
