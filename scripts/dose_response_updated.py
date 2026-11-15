# -*- coding: utf-8 -*-
"""
Created on Tue May 19 07:54:16 2026

@author: srahm
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ============================================================
# Mechanistic dose-response simulation
# ============================================================

np.random.seed(123)

# -----------------------------
# Baseline model parameters
# -----------------------------
params = {
    "lambda": 1.0,      # pathogen replication rate
    "K": 1e8,           # within-host carrying capacity
    "phi": 2e-5,        # immune-mediated killing
    "mu": 0.10,         # baseline pathogen loss
    "alpha": 1e-3,      # immune activation by pathogen
    "beta": 0.30,       # immune expansion rate
    "Imax": 1e5,        # immune carrying capacity
    "delta": 0.20,      # immune decay
    "I0": 4e4           # initial immune activity
}

# -----------------------------
# Simulation settings
# -----------------------------
doses = np.logspace(0, 8, 12)   # from 10^0 to 10^8
N_hosts = 500
T = 30
t_eval = np.linspace(0, T, 400)

# Host variability
sigma_eta = 0.4
sigma_phi = 0.25

# Mean survival/deposition fraction
eta_mean = 1e-3

# Establishment threshold
rho = 1e-2
Pc = rho * params["K"]

print(f"Establishment threshold Pc = {Pc:.3e}")
print(f"Carrying capacity K = {params['K']:.3e}")
print(f"Threshold fraction rho = {rho:.3e}")

# Representative doses for trajectory figure:
# low, middle, high
trajectory_doses = [doses[0], doses[len(doses)//2], doses[-1]]

# -----------------------------
# ODE system
# -----------------------------
def pathogen_immune_model(t, y, lam, K, phi, mu, alpha, beta, Imax, delta):
    P, I = y

    # numerical safety
    P = max(P, 0.0)
    I = max(I, 0.0)

    dPdt = lam * P * (1 - P / K) - phi * P * I - mu * P
    dIdt = alpha * P + beta * I * (1 - I / Imax) - delta * I

    return [dPdt, dIdt]


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

        # Lognormal survival/deposition variability
        eta_i = np.random.lognormal(
            mean=np.log(eta_mean) - 0.5 * sigma_eta**2,
            sigma=sigma_eta
        )

        # Lognormal immune-clearance variability
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

        P_traj = np.maximum(sol.y[0], 0.0)
        I_traj = np.maximum(sol.y[1], 0.0)

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

        # Save a few representative trajectories
        if any(np.isclose(D, td) for td in trajectory_doses) and i < 4:
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
        "threshold_Pc": Pc,
        "threshold_fraction_rho": rho,
        "N_hosts": N_hosts
    })

summary_df = pd.DataFrame(summary_rows)
host_df = pd.DataFrame(all_host_rows)

if trajectory_examples:
    trajectory_df = pd.concat(trajectory_examples, ignore_index=True)
else:
    trajectory_df = pd.DataFrame()

# ============================================================
# Save Excel file
# ============================================================
output_excel = "simulation_summary.xlsx"

with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="Dose_Response", index=False)
    host_df.to_excel(writer, sheet_name="Virtual_Hosts", index=False)
    trajectory_df.to_excel(writer, sheet_name="Example_Trajectories", index=False)

print(f"Saved Excel output: {output_excel}")

# ============================================================
# Figure 1: Emergent dose-response curve
# ============================================================
plt.figure(figsize=(6.5, 4.8))

plt.plot(
    summary_df["dose"],
    summary_df["infection_probability"],
    marker="o",
    linewidth=2.5,
    markersize=7,
    label="Mechanistic simulation"
)

plt.xscale("log")
plt.ylim(-0.05, 1.05)
plt.xlabel("Dose, $D$")
plt.ylabel(r"Estimated infection probability, $\widehat{P}_{\mathrm{inf}}(D)$")
plt.title("Emergent dose--response curve from mechanistic simulations")
plt.legend()
plt.tight_layout()
plt.savefig("dose_response_curve.png", dpi=300)
plt.close()

print("Saved figure: dose_response_curve.png")

# ============================================================
# Figure 2: Representative within-host trajectories
# ============================================================
dose_colors = {
    trajectory_doses[0]: "blue",
    trajectory_doses[1]: "orange",
    trajectory_doses[2]: "red"
}

plt.figure(figsize=(7.0, 5.0))

for D in trajectory_doses:
    subset = trajectory_df[np.isclose(trajectory_df["dose"], D)]

    first_label = True
    for host_id in subset["host_id"].unique():
        s = subset[subset["host_id"] == host_id]

        plt.plot(
            s["time"],
            s["pathogen_burden"],
            color=dose_colors[D],
            alpha=0.55,
            linewidth=1.5,
            label=f"$D={D:.1e}$" if first_label else None
        )
        first_label = False

plt.axhline(
    Pc,
    linestyle="--",
    linewidth=2,
    color="black",
    label=r"Operational threshold $P_c$"
)

plt.yscale("log")
plt.xlabel("Time")
plt.ylabel(r"Pathogen burden, $P(t)$")
plt.title("Representative within-host trajectories")
plt.legend()
plt.tight_layout()
plt.savefig("trajectory_examples.png", dpi=300)
plt.close()

print("Saved figure: trajectory_examples.png")

# ============================================================
# Print summary
# ============================================================
print("\nDose-response summary:")
print(summary_df[["dose", "infection_probability", "mean_max_pathogen_burden"]])