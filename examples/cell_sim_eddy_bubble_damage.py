import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from of_compartments import CompartmentModel, utils
from of_compartments.cellgrowth.damage_models import (
    bubbleDamageModel_cherry1992, bubbleDamageModel_walls2017,
    eddyDamageModel_lakhotia_papoutsakis)
from of_compartments.cellgrowth.growth_models import \
    growthModel_xing_simplified as growthModel
from of_compartments.cellgrowth.growth_models import \
    params_xing_simplified as params

plt.rcParams.update({"font.size": 14})


def plot_comparison(
    baseline_sol,
    other_sols,
    titles,
    labels,
    outfile,
):
    _, axs = plt.subplots(
        1, len(other_sols), layout="tight", figsize=(6 * len(other_sols), 5)
    )

    cell_counts = baseline_sol.z[:, 5, :] * baseline_sol.v
    vcd_baseline = np.sum(cell_counts, axis=1) / (1000 * np.sum(baseline_sol.v))
    vcd_baseline /= 1e6
    time = np.max(baseline_sol.t)
    for j in range(len(other_sols)):
        axs[j].plot(baseline_sol.t, vcd_baseline, "k:", label="no damage")
        for i in range(len(other_sols[j])):
            cell_counts = other_sols[j][i].z[:, 5, :] * other_sols[j][i].v
            vcd = np.sum(cell_counts, axis=1) / (1000 * np.sum(other_sols[j][i].v))
            vcd /= 1e6
            axs[j].plot(other_sols[j][i].t, vcd, label=labels[i])
            axs[j].set_ylabel("Viable Cell Density ($10^6$ cells / mL)")
            axs[j].set_title(titles[j], pad=20)

    for ax in axs.reshape(-1):
        ax.grid(axis="y")
        ax.set_xticks(np.arange(0, time + 1, 96))
        ax.set_xlim([0, time])
        ax.set_xlabel("Time (hr)")
        ax.legend()

    plt.savefig(outfile, format="svg")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a bioreactor sim.")
    parser.add_argument(
        "--case_names",
        "-c",
        nargs="+",
        help="names of the cases to use",
    )
    parser.add_argument(
        "--compartment_data_dir",
        "-d",
        default="compartment_data",
        help="directory containing the compartment data",
    )
    parser.add_argument(
        "--labels",
        "-l",
        nargs="+",
        help="labels of the cases to use",
    )
    parser.add_argument(
        "--time",
        "-t",
        default=244,
        help="time to run in hours",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        default=".",
        help="directory in which to write the output image",
    )

    args = parser.parse_args()
    case_names = args.case_names
    compartment_data_dir = args.compartment_data_dir
    case_labels = args.labels
    time = args.time

    X0 = 0.2e6  # cells / mL
    initial_concentrations = [100, 10, 0, 0, params["DO_eq"] / 5, X0 * 1000]

    eddy_sols = []
    bubble_walls_sols = []
    bubble_cherry_sols = []
    for case_name in case_names:
        case_dir = os.path.join(compartment_data_dir, case_name)
        compartments = utils.read_compartment_config(
            os.path.join(case_dir, "compartment_config")
        )
        cm_eddy = CompartmentModel(
            case_dir,
            damage_models=[
                eddyDamageModel_lakhotia_papoutsakis(
                    k_c=30, E_0=0.41, B=3e-4, nu=0.0007 / 993
                ),
            ],
        )
        sol_eddy = cm_eddy.run_sim(
            time=time,
            growth_model=growthModel,
            params=params,
            initial_concentrations=initial_concentrations,
        )
        eddy_sols.append(sol_eddy)

        cm_bubble_walls = CompartmentModel(
            case_dir,
            damage_models=[
                bubbleDamageModel_walls2017(
                    bubble_radius=0.002, threshold="mid", death_fraction=0.4
                ),
            ],
        )
        sol_bubble_walls = cm_bubble_walls.run_sim(
            time=time,
            growth_model=growthModel,
            params=params,
            initial_concentrations=initial_concentrations,
        )
        bubble_walls_sols.append(sol_bubble_walls)

        cm_bubble_cherry = CompartmentModel(
            case_dir,
            damage_models=[bubbleDamageModel_cherry1992(bubble_radius=0.002)],
        )
        sol_bubble_cherry = cm_bubble_cherry.run_sim(
            time=time,
            growth_model=growthModel,
            params=params,
            initial_concentrations=initial_concentrations,
        )
        bubble_cherry_sols.append(sol_bubble_cherry)

        cm_baseline = CompartmentModel(
            case_dir,
            damage_models=[],
        )
        sol_baseline = cm_baseline.run_sim(
            time=time,
            growth_model=growthModel,
            params=params,
            initial_concentrations=initial_concentrations,
        )

    if case_labels is None:
        case_labels = [name.split("_")[1] for name in case_names]

    plot_comparison(
        sol_baseline,
        other_sols=[eddy_sols, bubble_walls_sols, bubble_cherry_sols],
        titles=["Eddy damage", "Walls 2017 bubble damage", "Cherry 1992 bubble damage"],
        labels=case_labels,
        outfile=os.path.join(args.output_dir, "eddy_bubble_damage.svg"),
    )
