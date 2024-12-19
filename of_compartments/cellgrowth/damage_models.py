from collections.abc import Callable

import numpy as np
import scipy.integrate as si

from .compartment_model import CompartmentData


def bubbleDamageModel_shell(
    bubble_radius: float, cell_radius: float
) -> Callable[[CompartmentData], np.ndarray]:
    """
    Calculates damage rate based on liquid volume fraction that is
    within a cell radius of a bubble.

    From:
    Farzan & Ierapetritou (2018) A Framework for the Development of
    Integrated and Computationally Feasible Models of Large-Scale
    Mammalian Cell Bioreactors, Processes.
    """

    def bubbleDamageModel(compartment_data: CompartmentData) -> np.ndarray:
        bubble_lifespan = (
            compartment_data.gas_holdup
            * compartment_data.volumes
            / compartment_data.top_gas_flux
        )
        # Get the fraction of total volume that is in shells of thickness
        # equal to the cell radius around bubbles
        gas_volume = compartment_data.gas_holdup * compartment_data.volumes
        bubble_volume = (4 * np.pi / 3) * bubble_radius**3
        num_bubbles = gas_volume / bubble_volume
        interaction_volume = (4 * np.pi / 3) * (
            (bubble_radius + cell_radius) ** 3 - bubble_radius**3
        )
        bubble_danger_volume_fraction = (
            num_bubbles * interaction_volume
        ) / compartment_data.volumes
        # k_d_bubble
        return np.log(1 - bubble_danger_volume_fraction) / bubble_lifespan

    return bubbleDamageModel


def bubbleDamageModel_cherry1992(
    bubble_radius: float,
    Psi: float = 0.2,
    film_cell_density_ratio: float = 0.6,
    film_thickness: float = 1.5e-5,
) -> Callable[[CompartmentData], np.ndarray]:
    """
    bubble_radius: radius of bubbles in m
    Psi: fraction of cells that die in a bubble bursting event (default 20%)
    film_cell_density_ratio: ratio of cell density in thin film around bubble to
    that in liquid bulk (default 0.6)
    film_thickness: thickness of film around bubble in meters (default 1.5e-5 m = 15 um)

    Based on:
    Cherry & Hulle (1992) Cell Death in the Thin Films of Bursting Bubbles, Biotechnology Progress.
    """

    def bubbleDamageModel(compartment_data: CompartmentData) -> np.ndarray:
        top_compartments = np.array(compartment_data.ids)[compartment_data.is_top]
        top_indices = np.isin(compartment_data.ids, top_compartments)
        k = (
            Psi
            * film_cell_density_ratio
            * 3
            * film_thickness
            * compartment_data.top_gas_flux
            / (bubble_radius * compartment_data.volumes)
        )
        k *= 3600  # convert to h^-1

        k_d = np.zeros(len(compartment_data.ids))
        k_d[top_indices] = k[top_indices]
        return k_d

    return bubbleDamageModel


def bubbleDamageModel_walls2017(
    bubble_radius: float, threshold: str, death_fraction: float = 1
) -> Callable[[CompartmentData], np.ndarray]:
    """
    bubble_radius: radius of bubbles in m

    threshold: can be 'low' (10^6 W / m^3) or 'mid' (10^7 W/m^3)

    death_fraction: fraction of cells killed on exposure to threshold shear

    Based on data from Walls, McRae & Natarajan et al. (2017) Quantifying the
    Potential for Bursting Bubbles to Damage Suspended Cells, Scientific Reports.

    """
    # convert to um
    R = bubble_radius * 1e6

    if threshold == "low":
        m = -1.550
        b = 7.982
    elif threshold == "mid":
        m = -1.848
        b = 8.261
    else:
        raise Exception(f"invalid threshold value {threshold}")

    def bubbleDamageModel(compartment_data: CompartmentData) -> np.ndarray:
        top_compartments = np.array(compartment_data.ids)[compartment_data.is_top]
        top_indices = np.isin(compartment_data.ids, top_compartments)
        coef = np.exp(b) * R**m
        V_b = 4 / 3 * np.pi * (R * 1e-6) ** 3  # m^3
        V_deadly = V_b * coef
        num_bubbles = compartment_data.top_gas_flux / V_b
        V_deadly_total = V_deadly * num_bubbles
        q = 3600 * V_deadly_total / compartment_data.volumes

        k_d = np.zeros(len(compartment_data.ids))
        k_d[top_indices] = q[top_indices]
        return k_d * death_fraction

    return bubbleDamageModel


def shearDamageModel_constant(rate: float) -> Callable[[CompartmentData], np.ndarray]:
    def shearDamageModel(compartment_data: CompartmentData) -> np.ndarray:
        return rate * compartment_data.high_shear_fraction

    return shearDamageModel


def eddyDamageModel_lakhotia_papoutsakis(
    B, k_c, E_0, nu
) -> Callable[[CompartmentData], np.ndarray]:
    """
    Based on:
    Lakhotia & Papoutsakis (1992) Agitation Induced Cell Injury in Microcarrier Cultures.
    Protective Effect of Viscosity Is Agitation Intensity Dependent: Experiments and Modeling,
    Biotechnology and Bioengineering.

    B in s/cm^2
    k_c in cm^-1
    E_0 in cm^2/s^2
    """

    def eddyDamageModel(compartment_data: CompartmentData) -> np.ndarray:
        # convert from m^2/s^3 to cm^2/s^3
        eps = compartment_data.epsilon * 1e4

        I, _ = si.quad_vec(
            lambda k: 1.7
            * eps ** (2 / 3)
            * k ** (-5 / 3)
            * np.exp(-2.55 * nu * eps ** (-1 / 3) * k ** (4 / 3)),
            k_c,
            np.inf,
        )
        k_d = np.zeros(len(eps))
        k_d[I > E_0] = B * I[I > E_0]
        return k_d

    return eddyDamageModel
