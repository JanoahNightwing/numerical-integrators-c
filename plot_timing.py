import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load(path):
    rows = {"euler": [], "rk4": []}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            rows[r["method"]].append((float(r["h"]), float(r["error"]), float(r["seconds"])))
    return {k: np.array(v) for k, v in rows.items()}

o0, o2 = load("timing_O0.csv"), load("timing_O2.csv")
fig, ax = plt.subplots(1, 2, figsize=(11, 4))

# Left: error vs runtime at -O2, stable and above the rounding floor only.
for name, key in (("Euler", "euler"), ("RK4", "rk4")):
    h, err, sec = o2[key][:, 0], o2[key][:, 1], o2[key][:, 2]
    keep = (h < 0.1) & (err > 1e-14)
    ax[0].loglog(sec[keep], err[keep], "o-", label=name)
ax[0].set_xlabel("time per solve (s)")
ax[0].set_ylabel("error at t = 10")
ax[0].set_title("Error vs runtime (-O2)")
ax[0].grid(True, which="both")
ax[0].legend()

# Right: median nanoseconds per step, for runs with at least 1000 steps.
labels, vals = [], []
print("median ns per step (runs with >= 1000 steps):")
for opt, d in (("-O0", o0), ("-O2", o2)):
    for key in ("euler", "rk4"):
        n = np.rint(10.0 / d[key][:, 0])
        ns = d[key][:, 2] / n * 1e9
        med = float(np.median(ns[n >= 1000]))
        print(f"  {opt} {key}: {med:.1f}")
        labels.append(f"{opt}\n{key}")
        vals.append(med)
ax[1].bar(range(len(vals)), vals)
ax[1].set_xticks(range(len(vals)))
ax[1].set_xticklabels(labels)
ax[1].set_ylabel("ns per step")
ax[1].set_title("Cost per step")

plt.tight_layout()
plt.savefig("timing.png", dpi=150)