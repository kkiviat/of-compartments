import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
from of_compartments import CompartmentModel, utils
from of_compartments.cellgrowth.growth_models import \
    growthModel_xing_simplified as growthModel
from of_compartments.cellgrowth.growth_models import \
    params_xing_simplified as params


def plot_results(sols, labels, outfile):
    fig, axs = plt.subplots(2, 3, layout="tight", figsize=(8, 4))
    time = np.max(sols[0].t)
    for i, sol in enumerate(sols):
        cell_counts = sol.z[:, 5, :] * sol.v
        vcd = np.sum(cell_counts, axis=1) / (1000 * np.sum(sol.v))
        vcd /= 1e6
        axs[0, 0].plot(sol.t, vcd, label=labels[i])
        axs[0, 0].set_title("Viable Cell Density ($10^6$ cells / mL)", pad=20)

        counts = 1000 * sol.z[:, 0, :] * sol.v
        glucose = np.sum(counts, axis=1) / (1000 * np.sum(sol.v))
        axs[0, 1].plot(sol.t, glucose, label=labels[i])
        axs[0, 1].set_title("Glucose (mM)", pad=20)

        counts = 1000 * sol.z[:, 1, :] * sol.v
        glutamine = np.sum(counts, axis=1) / (1000 * np.sum(sol.v))
        axs[0, 2].plot(sol.t, glutamine, label=labels[i])
        axs[0, 2].set_title("Glutamine (mM)", pad=20)

        counts = 1000 * sol.z[:, 2, :] * sol.v
        lactate = np.sum(counts, axis=1) / (1000 * np.sum(sol.v))
        axs[1, 0].plot(sol.t, lactate, label=labels[i])
        axs[1, 0].set_title("Lactate (mM)", pad=20)

        counts = 1000 * sol.z[:, 3, :] * sol.v
        ammonia = np.sum(counts, axis=1) / (1000 * np.sum(sol.v))
        axs[1, 1].plot(sol.t, ammonia, label=labels[i])
        axs[1, 1].set_title("Ammonia (mM)", pad=20)

        counts = 1000 * sol.z[:, 4, :] * sol.v
        DO = np.sum(counts, axis=1) / (1000 * np.sum(sol.v))
        DO *= 100 / 1.07
        axs[1, 2].plot(sol.t, DO, label=labels[i])
        axs[1, 2].set_title("DO (%)", pad=20)
        axs[1, 2].set_ylim([0, 100])

    for ax in axs.reshape(-1):
        ax.grid(axis="y")
        ax.set_xticks(np.arange(0, time + 1, 96))
        ax.set_xlim([0, time])
        ax.set_xlabel("Time (hr)")

    if len(sols) > 1:
        handles, labels = axs[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper right")

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
        "--time",
        "-t",
        default=244,
        help="time to run in hours",
        required=True,
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

    X0 = 0.2e6  # cells / mL
    initial_concentrations = [100, 10, 0, 0, 1.07 / 5, X0 * 1000]

    sols = []
    for case_name in case_names:
        case_dir = os.path.join(compartment_data_dir, case_name)
        compartments = utils.read_compartment_config(
            os.path.join(case_dir, "compartment_config")
        )

        cm = CompartmentModel(
            case_dir,
            damage_models=[],
        )
        sol = cm.run_sim(
            time=args.time,
            growth_model=growthModel,
            params=params,
            initial_concentrations=initial_concentrations,
        )
        sols.append(sol)

    plot_results(sols, case_names, os.path.join(args.output_dir, "cell_sim_batch.svg"))
