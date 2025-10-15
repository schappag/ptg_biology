import numpy as np
import matplotlib.pyplot as plt

# Parameters
pD = [5e-4, 1e-3]   # [rate, sensitivity]
y_pat3 = 900        # baseline PTH
y21 = 1.0           # initial state value

# Original formulation
def dydt21_original(y3, y21=y21, pD=pD, y_pat3=y_pat3):
    return pD[0] * (1 + pD[1] * (y3 - y_pat3) * y21)

# Stabilized formulation with tanh
def dydt21_tanh(y3, y21=y21*50, pD=pD, y_pat3=y_pat3):
    scale_PTH = 65
    return pD[0] * (1 + pD[1] * np.tanh((y3 - y_pat3) / scale_PTH) * y21)


    dydt[21] = pD[0] * (1 + pD[1] * (y[3] - y_pat[3]) * y[21])
    dydt[22] = pDD[0] * (1 + pDD[1] * (D-D_pat) -  y[22])
    threshold = 1e-6


# Sweep values of y[3]
y3_values = np.linspace(0, 2000, 500)
orig_vals = [dydt21_original(y3) for y3 in y3_values]
tanh_vals = [dydt21_tanh(y3) for y3 in y3_values]

# Plot
plt.figure(figsize=(8,5))
plt.plot(y3_values, orig_vals, label="Original")
plt.plot(y3_values, tanh_vals, label="Stabilized (tanh)", linestyle="--")
plt.axvline(y_pat3, color="k", linestyle=":", label="baseline y[3]")
plt.xlabel("y[3] (plasma PTH)")
plt.ylabel("dydt[21]")
plt.legend()
plt.title("Comparison of original vs. tanh-stabilized dydt[21]")
plt.grid(True)
plt.show()
