import numpy as np
import matplotlib
matplotlib.use("Agg")  # no display needed; saves to a file
import matplotlib.pyplot as plt

d = np.loadtxt("convergence.csv", delimiter=",", skiprows=1)
h, e_err, r_err = d[:, 0], d[:, 1], d[:, 2]

def slope(h, err, mask):
    return np.polyfit(np.log(h[mask]), np.log(err[mask]), 1)[0]

# Fit only where each method is in its asymptotic regime.
e_slope = slope(h, e_err, h < 0.01)
r_slope = slope(h, r_err, (h < 0.2) & (r_err > 1e-11))  # above the rounding floor

plt.loglog(h, e_err, "o-", label=f"Euler (fitted slope {e_slope:.2f})")
plt.loglog(h, r_err, "s-", label=f"RK4 (fitted slope {r_slope:.2f})")
plt.xlabel("step size h")
plt.ylabel("error at t = 10")
plt.grid(True, which="both")
plt.legend()
plt.savefig("convergence.png", dpi=150)
print("euler slope:", e_slope)
print("rk4 slope:", r_slope)