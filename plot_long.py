import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.loadtxt("energy_long.csv", delimiter=",", skiprows=1)
t, rk, ve = d[:, 0], d[:, 1], d[:, 2]

plt.plot(t, rk - 1, label="RK4")
plt.plot(t, ve - 1, ".", markersize=2, label="velocity Verlet")
plt.xlabel("t")
plt.ylabel("E / E0 - 1")
plt.grid(True)
plt.legend()
plt.savefig("energy_long.png", dpi=150)