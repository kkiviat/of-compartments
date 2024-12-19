import os
from collections.abc import Callable
from dataclasses import dataclass
from typing import List

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


@dataclass
class CompartmentData:
    ids: List[str]
    gas_holdup: np.ndarray
    high_shear_fraction: np.ndarray
    kLa: np.ndarray
    volumes: np.ndarray
    top_gas_flux: np.ndarray
    epsilon: np.ndarray
    is_top: np.ndarray


@dataclass
class Case:
    volume: float
    F: np.ndarray
    k_d: np.ndarray
    compartment_data: CompartmentData


class CompartmentModel:
    def __init__(
        self,
        case_dir: str,
        damage_models: List[Callable[[CompartmentData], np.ndarray]],
    ):
        self.case_dir = case_dir
        self.damage_models = damage_models

        self.case = self._create_case()

    def run_sim(
        self,
        time: float,
        growth_model: Callable[[float, np.ndarray, dict, Case], np.ndarray],
        params: dict,
        initial_concentrations: List[float],
    ) -> None:
        """
        Run the cell growth simulation in the compartment model

        Args:
            time: The time to run the simulation for (hours)
            growthModel: Function describing the cell growth.
            params: dict specifying parameters needed by the growth model.
            initial_concentrations: array of initial values for cells and
                metabolites, in the same order as used in the growthModel.
        """
        num_compartments = len(self.compartment_ids)

        IV = np.zeros((len(initial_concentrations), num_compartments))
        for i in range(len(initial_concentrations)):
            IV[i, :] = num_compartments * [initial_concentrations[i]]

        IV = IV.reshape(
            IV.shape[0] * IV.shape[1],
        )

        sol = solve_ivp(
            CompartmentModel._update_model,
            [0, time],
            IV,
            args=(
                params,
                self.case,
                self.compartment_ids,
                growth_model,
            ),
            method="BDF",
            atol=1e-12,
            rtol=1e-8,
            # max_step=0.01,
        )

        self.t = sol.t
        self.x = sol.y
        z = sol.y.T
        sol.z = z.reshape(z.shape[0], len(initial_concentrations), num_compartments)
        sol.v = self.case.compartment_data.volumes
        return sol

    @staticmethod
    def _readCompartmentData(case_dir: str) -> CompartmentData:
        """Reads compartment data from the given case directory."""
        compartment_df = pd.read_csv(os.path.join(case_dir, "compartment_values.csv"))
        compartment_ids = compartment_df["compartment"].values
        gas_holdup = compartment_df["gas_holdup"].values
        kLa = compartment_df["kLa"].values
        epsilon = compartment_df["epsilon"].values
        volumes = compartment_df["volume"].values
        top_gas_flux = compartment_df["top_gas_flux"].values
        high_shear_fraction = compartment_df["high_tau_fraction"].values
        is_top = compartment_df["is_top"].values

        return CompartmentData(
            ids=compartment_ids,
            gas_holdup=gas_holdup,
            high_shear_fraction=high_shear_fraction,
            kLa=kLa,
            volumes=volumes,
            top_gas_flux=top_gas_flux,
            epsilon=epsilon,
            is_top=is_top,
        )

    def _create_case(self) -> Case:
        compartment_data = CompartmentModel._readCompartmentData(self.case_dir)
        self.compartment_ids = list(compartment_data.ids)
        k_d = np.zeros(len(self.compartment_ids))
        for model in self.damage_models:
            k_d += model(compartment_data)

        F = self._create_interface_values(self.case_dir, self.compartment_ids)
        return Case(
            F=F,
            compartment_data=compartment_data,
            k_d=k_d,
            volume=np.sum(compartment_data.volumes),
        )

    def _create_interface_values(
        self, case_dir: str, compartment_ids: List[str]
    ) -> np.ndarray:
        """Calculates flow matrix between compartments in m^3/h."""
        F = np.zeros((len(compartment_ids), len(compartment_ids)))
        interface_values = pd.read_csv(
            os.path.join(case_dir, "interface_values.csv"), index_col=None
        )
        for _, row in interface_values.iterrows():
            c1 = compartment_ids.index(row["compartment_src"])
            c2 = compartment_ids.index(row["compartment_dest"])
            F[c1, c2] = row["corrected_flow"] * 3600  # convert m^3/s -> m^3/h

        return F

    @staticmethod
    def _get_flux_array(
        concentrations: np.ndarray, F: np.ndarray, volumes: np.ndarray
    ) -> np.ndarray:
        """Converts concentrations to a numerical flux matrix between all pairs of compartments."""
        # matrix [c_i * F_ij]
        X = F * concentrations[:, np.newaxis]
        # num entering compartment i is sum of column i of X
        # num leaving compartment i is sum of row i of X
        return (np.sum(X, axis=0) - np.sum(X, axis=1)) / volumes

    @staticmethod
    def _update_model(
        t: float,
        z: np.ndarray,
        p: dict,
        case: Case,
        compartment_ids: List[str],
        grow_cells: Callable[[float, np.ndarray, dict, Case], np.ndarray],
    ) -> np.ndarray:
        """
        Solves given growth model in compartments, handling flux between compartments.

        Args:
            t: Time to solve growth model (hours)
            z: Array of concentrations of cells and metabolites in each compartment
                (rows correspond to cells / metabolites, columns to compartments.
            p: Dictionary of parameter values for growCells.
            case: Case data
            compartment_ids: list of compartment id strings
            growCells: Function describing the cell growth.
        """
        num_concentrations = int(len(z) / len(compartment_ids))
        z = z.reshape(num_concentrations, len(compartment_ids))
        growth_deltas = grow_cells(t, z, p, case)
        F = case.F
        volumes = case.compartment_data.volumes

        result = np.zeros(z.shape)
        for i in range(z.shape[0]):
            result[i, :] = growth_deltas[i, :] + CompartmentModel._get_flux_array(
                z[i, :], F, volumes
            )
        result = result.reshape(
            result.shape[0] * result.shape[1],
        )

        return result
