# -*- coding: utf-8 -*-
"""
Created on Sun May 17 21:00:54 2026

@author: srahm
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

# ============================================================
# Mechanistic dose-response simulation
# ============================================================

np.random.seed(123)

# -----------------------------
# Baseline model parameters
# -----------------------------
params = {
    "lambda": 1.0,      # pathogen replication rate
    "K": 1e8,           # carrying capacity
    "phi": 2e-5,        # immune-mediated killing
    "mu": 0.10,         # baseline pathogen loss
    "alpha": 1e-3,      # immune activation by pathogen
    "beta": 0.30,       # immune expansion rate
    "Imax": 1e5,        # immune carrying capacity
    "delta": 0.20,      # immune decay
    "I0": 4e4          # initial immune activity
}

# -----------------------------
# Simulation settings
# -----------------------------
#doses = np.array([10**k for k in range(0, 7)], dtype=float)
doses = np.logspace(0, 8, 10)
N_hosts = 500
T = 30
t_eval = np.linspace(0, T, 400)

# Host variability
sigma_eta = 0.4
sigma_phi = 0.25

# Mean survival fraction
eta_mean = 1e-3

# -----------------------------
# Mechanistic threshold Pc
# -----------------------------
lam = params["lambda"]
mu = params["mu"]
phi0 = params["phi"]
alpha = params["alpha"]
beta = params["beta"]
Imax = params["Imax"]
delta = params["delta"]

Ic = (lam - mu) / phi0
rho = 1e-2
Pc = rho * params["K"]
# =============================================================================
# Pc = (
#     delta * Ic
#     - beta * Ic * (1 - Ic / Imax)
# ) / alpha
# =============================================================================

if Ic >= Imax:
    raise ValueError(
        "Invalid threshold regime: Ic >= Imax. "
        "Increase phi or Imax so that Ic < Imax."
    )

if Pc <= 0 or Pc >= params["K"]:
    raise ValueError(
        f"Invalid Pc={Pc:.3e}. Choose parameters so 0 < Pc < K."
    )

print(f"Critical immune level Ic = {Ic:.3e}")
print(f"Mechanistic pathogen threshold Pc = {Pc:.3e}")

# -----------------------------
# ODE system
# -----------------------------
def pathogen_immune_model(t, y, lam, K, phi, mu, alpha, beta, Imax, delta):
    P, I = y

    P = max(P, 0)
    I = max(I, 0)

    dPdt = lam * P * (1 - P / K) - phi * P * I - mu * P
    dIdt = alpha * P + beta * I * (1 - I / Imax) - delta * I

    return [dPdt, dIdt]


# -----------------------------
# Exponential model for comparison
# -----------------------------
def exponential_model(D, r):
    return 1 - np.exp(-r * D)


# -----------------------------
# Storage
# -----------------------------
summary_rows = []
all_host_rows = []
trajectory_examples = []

# -----------------------------
# Main simulation
# -----------------------------
for D in doses:

    infected_count = 0
    max_burdens = []

    for i in range(N_hosts):

        # Lognormal variability
        eta_i = np.random.lognormal(
            mean=np.log(eta_mean) - 0.5 * sigma_eta**2,
            sigma=sigma_eta
        )

        phi_i = np.random.lognormal(
            mean=np.log(params["phi"]) - 0.5 * sigma_phi**2,
            sigma=sigma_phi
        )

        P0 = eta_i * D
        I0 = params["I0"]

        sol = solve_ivp(
            pathogen_immune_model,
            [0, T],
            [P0, I0],
            t_eval=t_eval,
            args=(
                params["lambda"],
                params["K"],
                phi_i,
                params["mu"],
                params["alpha"],
                params["beta"],
                params["Imax"],
                params["delta"],
            ),
            method="RK45",
            rtol=1e-6,
            atol=1e-9
        )

        P_traj = np.maximum(sol.y[0], 0)
        I_traj = np.maximum(sol.y[1], 0)

        Ymax = np.max(P_traj)
        infected = int(Ymax >= Pc)

        infected_count += infected
        max_burdens.append(Ymax)

        all_host_rows.append({
            "dose": D,
            "host_id": i + 1,
            "eta": eta_i,
            "phi": phi_i,
            "P0": P0,
            "max_pathogen_burden": Ymax,
            "infected": infected
        })

        # Save representative trajectories for selected doses
        if D in [1, 1e3, 1e6] and i < 4:
            trajectory_examples.append(pd.DataFrame({
                "time": sol.t,
                "pathogen_burden": P_traj,
                "immune_activity": I_traj,
                "dose": D,
                "host_id": i + 1
            }))

    Pinf = infected_count / N_hosts

    summary_rows.append({
        "dose": D,
        "infection_probability": Pinf,
        "mean_max_pathogen_burden": np.mean(max_burdens),
        "median_max_pathogen_burden": np.median(max_burdens),
        "establishment threshold $P_c$": Pc,
        "N_hosts": N_hosts
    })

summary_df = pd.DataFrame(summary_rows)
host_df = pd.DataFrame(all_host_rows)
trajectory_df = pd.concat(trajectory_examples, ignore_index=True)

# -----------------------------
# Fit exponential model
# -----------------------------
try:
    popt, _ = curve_fit(
        exponential_model,
        summary_df["dose"].values,
        summary_df["infection_probability"].values,
        bounds=(0, np.inf),
        maxfev=10000
    )
    r_fit = popt[0]
except Exception:
    r_fit = np.nan

summary_df["exponential_fit"] = exponential_model(summary_df["dose"], r_fit)

print(f"Fitted exponential r = {r_fit:.3e}")

# ============================================================
# Save Excel file
# ============================================================
output_excel = "simulation_summary.xlsx"

with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="Dose_Response", index=False)
    host_df.to_excel(writer, sheet_name="Virtual_Hosts", index=False)
    trajectory_df.to_excel(writer, sheet_name="Example_Trajectories", index=False)

print(f"Saved Excel output: {output_excel}")
plt.figure(figsize=(7.0, 5.0))
dose_colors = {
    1: "blue",
    1e3: "orange",
    1e6: "red"
}
color=dose_colors[D] # for color legend
# ============================================================
# Figure 1: Dose-response curve
# ============================================================
dose_grid = np.logspace(0, 6, 300)

plt.figure(figsize=(6.5, 4.8))
plt.plot(
    summary_df["dose"],
    summary_df["infection_probability"],
    marker="o",
    linewidth=2,
    label="Mechanistic simulation"
)

# =============================================================================
# if not np.isnan(r_fit):
#     plt.plot(
#         dose_grid,
#         exponential_model(dose_grid, r_fit),
#         linestyle="--",
#         linewidth=2,
#         label="Exponential fit"
#     )
# =============================================================================

plt.xscale("log")
plt.ylim(-0.05, 1.05)
plt.xlabel("Dose, $D$")
plt.ylabel(r"Estimated infection probability,$\widehat{P}_{\mathrm{inf}}(D)$")
plt.title("Emergent dose--response curve from mechanistic simulations")
plt.legend()
plt.tight_layout()
plt.savefig("dose_response_curve.png", dpi=300)
plt.close()

print("Saved figure: dose_response_curve.png")

# ============================================================
# Figure 2: Example trajectories
# ============================================================


for D in [1, 1e3, 1e6]:
    subset = trajectory_df[trajectory_df["dose"] == D]
    for host_id in subset["host_id"].unique():
        s = subset[subset["host_id"] == host_id]
        plt.plot(
            s["time"],
            s["pathogen_burden"],
            alpha=0.45,
            linewidth=1.2,
            label=f"D={D:.0e}" if host_id == subset["host_id"].unique()[0] else None
        )

plt.axhline(Pc, linestyle="--", linewidth=2, label=r"Threshold $P_c$")
plt.yscale("log")
plt.xlabel("Time")
plt.ylabel("Pathogen burden, $P(t)$")
plt.title("Representative within-host trajectories")
plt.legend()
plt.tight_layout()
plt.savefig("trajectory_examples.png", dpi=300)
plt.close()

print("Saved figure: trajectory_examples.png")