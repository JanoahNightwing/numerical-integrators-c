import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.loadtxt("energy.csv", delimiter=",", skiprows=1)
t, eu, rk, ve = d[:, 0], d[:, 1], d[:, 2], d[:, 3]

fig, ax = plt.subplots(1, 2, figsize=(11, 4))

ax[0].semilogy(t, eu)
ax[0].set_title("Forward Euler")
ax[0].set_xlabel("t")
ax[0].set_ylabel("E / E0")
ax[0].grid(True, which="both")

ax[1].plot(t, rk - 1, label="RK4")
ax[1].plot(t, ve - 1, label="velocity Verlet")
ax[1].set_title("RK4 vs velocity Verlet")
ax[1].set_xlabel("t")
ax[1].set_ylabel("E / E0 - 1")
ax[1].grid(True)
ax[1].legend()

plt.tight_layout()
plt.savefig("energy.png", dpi=150)