import re
import subprocess
from collections import namedtuple
from functools import partial

import numpy as np
from scipy.optimize import minimize


def _parse_field(field, filename):
    return subprocess.check_output(
        f"awk -F ':' '/^{field}/{{print $2}}' {filename}", shell=True
    ).strip()


Compartment = namedtuple(
    "Compartment",
    ["bottom_height", "top_height", "inner_radius", "outer_radius", "id"],
)


def read_compartment_config(filename):
    """
    Reads the given file to get the top heights and outer radii of compartments
    """
    compartments = []
    heights, radii = read_compartment_bounds(filename)
    heights = [0] + heights
    radii = [0] + radii
    for i in range(len(heights) - 1):
        h_bottom = heights[i]
        h_top = heights[i + 1]
        for j in range(len(radii) - 1):
            r_inner = radii[j]
            r_outer = radii[j + 1]
            compartments.append(
                Compartment(
                    top_height=h_top,
                    bottom_height=h_bottom,
                    inner_radius=r_inner,
                    outer_radius=r_outer,
                    id=f"h{i}r{j}",
                    # id=f"i{i}" if j == 0 else f"o{i}",
                )
            )

    return compartments


def read_compartment_bounds(filename):
    """
    Reads the given file to get the top heights and outer radii of compartments
    """
    heights = _parse_field("heights", filename)
    heights = list(map(float, heights.split(b" ")))

    radii = _parse_field("radii", filename)
    radii = list(map(float, radii.split(b" ")))

    return heights, radii


def get_top_compartments(compartments):
    """Returns the ids of the top compartments in the list"""
    max_level = max([get_height_index_from_compartment_id(c.id) for c in compartments])
    return [
        c
        for c in compartments
        if get_height_index_from_compartment_id(c.id) == max_level
    ]


def get_top_compartments_by_id(compartment_ids):
    """Returns the ids of the top compartments in the list"""
    max_level = max(
        [get_height_index_from_compartment_id(id) for id in compartment_ids]
    )
    return [
        id
        for id in compartment_ids
        if get_height_index_from_compartment_id(id) == max_level
    ]


def get_compartment_lower_upper_pairs(compartments):
    """Returns pairs of compartments (c1, c2) such that c1 is directly below c2"""
    compartment_ids = get_ids(compartments)
    return get_compartment_lower_upper_pairs_by_name(compartment_ids)


def get_compartment_lower_upper_pairs_by_name(compartment_ids):
    """Returns pairs of compartments (c1, c2) such that c1 is directly below c2"""
    pairs = []
    max_level = max(
        [get_height_index_from_compartment_id(id) for id in compartment_ids]
    )
    max_radius = max(
        [get_radius_index_from_compartment_id(id) for id in compartment_ids]
    )
    for h in range(max_level):
        for r in range(max_radius + 1):
            pairs.append((get_zone_id(h, r), get_zone_id(h + 1, r)))

    return pairs


def get_ids(compartments):
    return [c.id for c in compartments]


def get_zone_id(h, r):
    return f"h{h}r{r}"


def get_height_index_from_compartment_id(id):
    return int(re.match(r".*h(\d+)", id).groups()[0])


def get_radius_index_from_compartment_id(id):
    return int(re.match(r".*r(\d+)", id).groups()[0])


# match i in h<i>r<j>


def order_compartments(c1, c2):
    """
    Create a consistent ordering for pairs of compartments.

    Returns the lower one first, and if tied, the inner one
    """
    c1_height = get_height_index_from_compartment_id(c1)
    c2_height = get_height_index_from_compartment_id(c2)
    c1_radius = get_radius_index_from_compartment_id(c1)
    c2_radius = get_radius_index_from_compartment_id(c2)
    if c1_height < c2_height:
        return c1, c2
    elif c2_height < c1_height:
        return c2, c1

    if c1_radius < c2_radius:
        return c1, c2

    return c2, c1


def create_boundary_name(c1, c2):
    c1, c2 = order_compartments(c1, c2)
    return f"{c1}_{c2}"


def correct_flow(F, tol=1e-6):
    """
    Corrects the flow matrix F.

    Adds small adjustments to minimize the net flow in each compartment.
    """
    F_flat = F.ravel()
    nonzero_idx = np.where(F_flat != 0)
    F_nonzero = F_flat[nonzero_idx]

    def objective(eps):
        return np.sum(np.square(eps))

    def flow_equality_constraint(eps):
        eps_full = np.zeros(len(F_flat))
        np.put(eps_full, nonzero_idx, eps)
        eps = np.reshape(eps_full, (F.shape))
        M = F - eps
        diffs = M.sum(axis=1) - M.sum(axis=0)
        return np.sum(np.abs(diffs))

    # See https://stackoverflow.com/a/6805307
    def correction_bounds_constraint(i, eps):
        # want abs(eps_ij) <= abs(F_ij), i.e.,
        # abs(F_ij) - abs(eps_ij) >= 0
        return np.abs(F_nonzero[i]) - np.abs(eps[i])

    constraints = [
        {"type": "eq", "fun": flow_equality_constraint},
    ]

    for i in range(len(F_nonzero)):
        constraints.append(
            {"type": "ineq", "fun": partial(correction_bounds_constraint, i)}
        )

    eps_init = np.zeros(len(F_nonzero))
    sol = minimize(
        objective,
        eps_init,
        constraints=constraints,
        method="SLSQP",
        options={"ftol": tol, "maxiter": 1000},
    )
    eps = np.zeros(F.shape)
    np.put(eps, nonzero_idx, sol.x)
    eps = eps.reshape(F.shape)
    corrected_F = F - eps
    if objective(sol.x) > 1 or np.isnan(objective(sol.x)):
        print(f"objective too high. Trying again with increased tol {tol*10}")
        corrected_F = correct_flow(F, tol * 10)
    return corrected_F
