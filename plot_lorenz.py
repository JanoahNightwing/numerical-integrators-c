import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.loadtxt("lorenz.csv", delimiter=",", skiprows=1)
t, sb, sc = d[:, 0], d[:, 1], d[:, 2]

# Rough growth rate: fit the stretch where the separation is growing cleanly.
mask = (sb > 1e-8) & (sb < 1.0)
rate = np.polyfit(t[mask], np.log(sb[mask]), 1)[0]

plt.semilogy(t, sb, label=f"1e-10 perturbation (fitted growth rate {rate:.2f})")
plt.semilogy(t, sc, label="step size h vs h/2")
plt.xlabel("t")
plt.ylabel("distance between trajectories")
plt.grid(True, which="both")
plt.legend()
plt.savefig("lorenz.png", dpi=150)
print("fitted growth rate:", rate)