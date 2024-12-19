"""Visualizes 2D slice of compartments with average field values"""
import argparse
import os

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle
from of_compartments import utils

matplotlib.rcParams["font.size"] = 22


def main():
    args = _parse_args()
    field = args.field

    case_dir = args.case_dir
    compartment_config = args.compartment_config
    if not compartment_config:
        compartment_config = os.path.join(case_dir, "compartment_config")

    df = pd.read_csv(compartment_config)

    if field == "shear":
        field = "tau"
        title = "average shear stress (Pa)"
    elif field == "high_shear_fraction":
        field = "high_tau_fraction"
        title = "fraction above shear threshold"
    elif field == "gas_holdup":
        title = "gas holdup (volume fraction)"
    elif field == "kLa":
        title = "kLa (1/h)"
    elif field == "epsilon":
        title = "epsilon (W/m^3)"
    elif field == "high_shear_fraction":
        title = "fraction above shear threshold"

    values = dict(zip(df["compartment"], df[field]))

    heights, radii = utils.read_compartment_bounds(compartment_config)
    heights = [0] + heights

    draw_compartments(
        heights,
        radii,
        values,
        title,
        show_text=args.show_text,
        scale_max=args.scale_max,
        show_axes=args.show_axes,
    )


def draw_compartments(
    heights,
    radii,
    values,
    title,
    show_text=False,
    scale_max=None,
    show_axes=False,
):
    fig, ax = plt.subplots()
    padding = 1.25
    min_val = min(values.values())
    max_val = max(values.values()) if not scale_max else scale_max
    norm = matplotlib.colors.Normalize(vmin=min_val, vmax=max_val)
    cmap = matplotlib.colormaps["turbo"]
    for i in range(len(heights) - 1):
        for r in reversed(range(len(radii))):
            ax.add_patch(
                Rectangle(
                    (-radii[r], heights[i]),
                    2 * radii[r],
                    heights[i + 1] - heights[i],
                    edgecolor="black",
                    facecolor=cmap(norm(values[utils.get_zone_id(h=i, r=r)])),
                )
            )
            if show_text:
                x_i = -radii[r]
                y_i = (heights[i + 1] - heights[i]) * r / len(radii) + heights[i]
                plt.text(x_i, y_i, values[f"h{i}r{r}"], color="white")

    ax.set_xlim([-radii[-1] * padding, radii[-1] * padding])
    ax.set_ylim([0, max(heights) + radii[-1] * (padding - 1)])
    ax.set_aspect("equal")
    ax.set_title(title)
    fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
    if not show_axes:
        plt.axis("off")

    plt.ioff()
    plt.show()


def _parse_args():
    """Parse the commandline arguments"""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case_dir",
        "-c",
        help="case directory",
        required=True,
    )
    parser.add_argument(
        "--compartment_config",
        help="path to compartment config",
    )
    parser.add_argument(
        "--scale_max",
        help="maximum value of cmap",
        type=float,
    )
    parser.add_argument(
        "--field",
        help="which field to graph",
        choices=["shear", "gas_holdup", "kLa", "high_shear_fraction", "epsilon"],
        default="shear",
    )
    parser.add_argument(
        "--show_text",
        help="whether to plot text",
        action="store_true",
    )
    parser.add_argument(
        "--show_axes",
        help="whether to show axes with spatial dimensions",
        action="store_true",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
