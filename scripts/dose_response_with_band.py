# -*- coding: utf-8 -*-
"""
Created on Tue May 19 08:58:57 2026

@author: srahm
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ============================================================
# Mechanistic dose-response simulation with uncertainty band
# ============================================================

np.random.seed(123)

# -----------------------------
# Baseline model parameters
# -----------------------------
params = {
    "lambda": 1.0,
    "K": 1e8,
    "phi": 2e-5,
    "mu": 0.10,
    "alpha": 1e-3,
    "beta": 0.30,
    "Imax": 1e5,
    "delta": 0.20,
    "I0": 4e4
}

# -----------------------------
# Simulation settings
# -----------------------------
doses = np.logspace(0, 8, 12)
M_runs = 100
N_hosts = 300
T = 30
t_eval = np.linspace(0, T, 400)

sigma_eta = 0.4
sigma_phi = 0.25
eta_mean = 1e-3

rho = 1e-2
Pc = rho * params["K"]

trajectory_doses = [doses[0], doses[len(doses)//2], doses[-1]]

print(f"Establishment threshold Pc = {Pc:.3e}")
print(f"Dose range: {doses[0]:.1e} to {doses[-1]:.1e}")

# -----------------------------
# ODE model
# -----------------------------
def pathogen_immune_model(t, y, lam, K, phi, mu, alpha, beta, Imax, delta):
    P, I = y
    P = max(P, 0.0)
    I = max(I, 0.0)

    dPdt = lam * P * (1 - P / K) - phi * P * I - mu * P
    dIdt = alpha * P + beta * I * (1 - I / Imax) - delta * I

    return [dPdt, dIdt]

# -----------------------------
# Storage
# -----------------------------
dose_response_runs = []
all_host_rows = []
trajectory_examples = []

# -----------------------------
# Monte Carlo simulation
# -----------------------------
for m in range(M_runs):

    run_probs = []

    for D in doses:

        infected_count = 0

        for i in range(N_hosts):

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

            P_traj = np.maximum(sol.y[0], 0.0)
            I_traj = np.maximum(sol.y[1], 0.0)

            Ymax = np.max(P_traj)
            infected = int(Ymax >= Pc)

            infected_count += infected

            all_host_rows.append({
                "run": m + 1,
                "dose": D,
                "host_id": i + 1,
                "eta": eta_i,
                "phi": phi_i,
                "P0": P0,
                "max_pathogen_burden": Ymax,
                "infected": infected
            })

            # Save example trajectories only from first run
            if (
                m == 0
                and any(np.isclose(D, td) for td in trajectory_doses)
                and i < 4
            ):
                trajectory_examples.append(pd.DataFrame({
                    "time": sol.t,
                    "pathogen_burden": P_traj,
                    "immune_activity": I_traj,
                    "dose": D,
                    "host_id": i + 1
                }))

        run_probs.append(infected_count / N_hosts)

    dose_response_runs.append(run_probs)

dose_response_runs = np.array(dose_response_runs)

# -----------------------------
# Summary statistics
# -----------------------------
summary_df = pd.DataFrame({
    "dose": doses,
    "infection_probability_mean": dose_response_runs.mean(axis=0),
    "infection_probability_median": np.median(dose_response_runs, axis=0),
    "infection_probability_lower_95": np.percentile(dose_response_runs, 2.5, axis=0),
    "infection_probability_upper_95": np.percentile(dose_response_runs, 97.5, axis=0),
    "threshold_Pc": Pc,
    "threshold_fraction_rho": rho,
    "M_runs": M_runs,
    "N_hosts_per_run": N_hosts
})

runs_df = pd.DataFrame(
    dose_response_runs,
    columns=[f"Dose_{D:.3e}" for D in doses]
)
runs_df.insert(0, "run", np.arange(1, M_runs + 1))

host_df = pd.DataFrame(all_host_rows)

trajectory_df = (
    pd.concat(trajectory_examples, ignore_index=True)
    if trajectory_examples else pd.DataFrame()
)

# -----------------------------
# Save Excel file
# -----------------------------
output_excel = "simulation_summary_with_band.xlsx"

with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="Dose_Response_Summary", index=False)
    runs_df.to_excel(writer, sheet_name="Monte_Carlo_Runs", index=False)
    host_df.to_excel(writer, sheet_name="Virtual_Hosts", index=False)
    trajectory_df.to_excel(writer, sheet_name="Example_Trajectories", index=False)

print(f"Saved Excel output: {output_excel}")

# ============================================================
# Figure 1: Dose-response curve with uncertainty band
# ============================================================
plt.figure(figsize=(6.5, 4.8))

plt.fill_between(
    summary_df["dose"],
    summary_df["infection_probability_lower_95"],
    summary_df["infection_probability_upper_95"],
    alpha=0.25,
    label="95% simulation interval"
)

plt.plot(
    summary_df["dose"],
    summary_df["infection_probability_median"],
    marker="o",
    linewidth=2.5,
    markersize=7,
    label="Median mechanistic simulation"
)

plt.xscale("log")
plt.ylim(-0.05, 1.05)
plt.xlabel("Dose, $D$")
plt.ylabel(r"Estimated infection probability, $\widehat{P}_{\mathrm{inf}}(D)$")
plt.title("Emergent dose--response curve from mechanistic simulations")
plt.legend()
plt.tight_layout()
plt.savefig("dose_response_curve_with_band.png", dpi=300)
plt.close()

print("Saved figure: dose_response_curve_with_band.png")

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
    label=r"Establishment threshold $P_c$"
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

# -----------------------------
# Print summary
# -----------------------------
print("\nDose-response summary:")
print(summary_df[
    [
        "dose",
        "infection_probability_median",
        "infection_probability_lower_95",
        "infection_probability_upper_95"
    ]
])