import matplotlib.pyplot as plt

def plot_iPTH(sn, unit_iPTH="pg/mL", abs_rel="abs"):
    if abs_rel == "abs":
        plt.figure(figsize=(6, 4))
        plt.plot(sn.t / 24, sn.y[3, :] / 3 * 9.43)
        plt.xlabel("Time [days]")
        plt.ylabel(f"iPTH [{unit_iPTH}]")
        plt.title("iPTH")
    else:
        plt.figure(figsize=(6, 4))
        plt.plot(sn.t / 24, sn.y[3, :] / sn.y[3, 0])
        plt.xlabel("Time [days]")
        plt.ylabel("Relative iPTH")
        plt.title("iPTH")

def plot_states(sn, indices=range(0,4)):
    for i in indices:
        plt.figure(figsize=(6, 4))
        plt.plot(sn.t / 24, sn.y[i, :])
        plt.title(f'{i}')
