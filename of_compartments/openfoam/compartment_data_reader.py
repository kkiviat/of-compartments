"""Reads data from openfoam function object output"""
import os
import subprocess

import numpy as np
import pandas as pd
from of_compartments import utils

from . import COMPARTMENT_AVERAGE_FIELDS

_FLOW_FORMAT = "flux_data_{phase}.csv"
_VOL_AVG_FORMAT = "{case_dir}/postProcessing/volAvg_{zone}/{t}/volFieldValue.dat"
_AREA_FORMAT = (
    "{case_dir}/postProcessing/area_{zone1}_{zone2}/{t}/surfaceFieldValue.dat"
)
# not in postProcessing
_TAU_FORMAT = "{case_dir}/{t}/tau_s_cellZone-{zone}_zone"
_CELL_VOLUME_FORMAT = "{case_dir}/{t}/V_cellZone-{zone}_zone"


def _read_area_from_file(file_name):
    # Area is third line in file
    file_name = file_name.replace("(", "\\(").replace(")", "\\)")
    area = float(
        subprocess.check_output(f"sed '3q;d' {file_name}", shell=True)
        .strip()
        .split(b" ")[-1]
    )
    return area


def _read_vol_avg_field(zone, time, case_dir, field):
    file_name = _VOL_AVG_FORMAT.format(zone=zone, t=time, case_dir=case_dir)
    vol_avg = pd.read_csv(file_name, sep="\t", skiprows=3)
    vol_avg.columns = ["Time"] + COMPARTMENT_AVERAGE_FIELDS
    return float(vol_avg.iloc[0][field])


def read_gas_holdup(zone, time, case_dir):
    return _read_vol_avg_field(zone, time, case_dir, "alphaMean.air")


def read_kLa(zone, time, case_dir):
    return _read_vol_avg_field(zone, time, case_dir, "kLa")


def read_tau_average(zone, time, case_dir):
    return _read_vol_avg_field(zone, time, case_dir, "tau_s")


def read_epsilon(zone, time, case_dir):
    return _read_vol_avg_field(zone, time, case_dir, "epsilonMean.water")


def read_top_gas_flux(zones, case_dir):
    "Returns volumetric flow of gas through top of each compartment"

    file_name = os.path.join(case_dir, _FLOW_FORMAT.format(phase="air"))
    BOUNDARY_ZONE_FORMAT = "boundary_{zone1}_{zone2}_zone"
    BOUNDARY_ZONE_TOP_FORMAT = "boundary_{zone}_top_zone"
    top_gas_flux = np.zeros(len(zones))

    try:
        flow_rate = pd.read_csv(file_name, index_col=3)
        for zone1, zone2 in utils.get_compartment_lower_upper_pairs_by_name(zones):
            name = BOUNDARY_ZONE_FORMAT.format(zone1=zone1, zone2=zone2)
            top_gas_flux[zones.index(zone1)] = flow_rate.loc[name]["flux_pos"]

        for top_zone in utils.get_top_compartments_by_id(zones):
            name = BOUNDARY_ZONE_TOP_FORMAT.format(zone=top_zone)
            top_gas_flux[zones.index(top_zone)] = flow_rate.loc[name]["flux_pos"]

    except FileNotFoundError:
        print("no flow found")

    return top_gas_flux


def read_flow(zones, case_dir, phase):
    "Returns volumetric flow"
    file_name = os.path.join(case_dir, _FLOW_FORMAT.format(phase=phase))
    BOUNDARY_ZONE_FORMAT = "boundary_{zone1}_{zone2}_zone"
    flow_matrix = np.zeros((len(zones), len(zones)))

    try:
        flow_rate = pd.read_csv(file_name, index_col=3)

        for i in range(len(zones)):
            for j in range(len(zones)):
                if i == j:
                    continue

                forward_name = BOUNDARY_ZONE_FORMAT.format(
                    zone1=zones[i], zone2=zones[j]
                )
                backward_name = BOUNDARY_ZONE_FORMAT.format(
                    zone1=zones[j], zone2=zones[i]
                )
                if forward_name in flow_rate.index:
                    flow_matrix[i, j] = flow_rate.loc[forward_name]["flux_pos"]
                    flow_matrix[j, i] = abs(flow_rate.loc[forward_name]["flux_neg"])
                elif backward_name in flow_rate.index:
                    flow_matrix[i, j] = abs(flow_rate.loc[backward_name]["flux_neg"])
                    flow_matrix[j, i] = flow_rate.loc[backward_name]["flux_pos"]
    except FileNotFoundError:
        print("no flow found")

    return flow_matrix


def read_volume(zone, time, case_dir):
    """Returns volume of given zone (m^3)"""
    file_name = _VOL_AVG_FORMAT.format(zone=zone, t=time, case_dir=case_dir)
    # Volume is third line in file
    volume = float(
        subprocess.check_output(f"sed '3q;d' {file_name}", shell=True)
        .strip()
        .split(b" ")[-1]
    )
    return volume


def read_cell_volumes(zone, time, case_dir):
    """Returns list of volumes of all cells in zone (m^3)"""
    file_name = _CELL_VOLUME_FORMAT.format(zone=zone, t=time, case_dir=case_dir)
    volumes = pd.read_fwf(file_name, skiprows=18).iloc[:-2]
    volumes = list(map(float, volumes["("]))
    return volumes


def read_tau_threshold_fraction(zone, time, threshold, case_dir):
    """Returns the volume fraction where shear stress is above threshold"""
    file_name = _TAU_FORMAT.format(t=time, zone=zone, case_dir=case_dir)
    shear = pd.read_csv(file_name, skiprows=18).iloc[:-2]
    shear = list(map(float, shear["("]))
    cells_above_threshold = list(map(lambda x: 1 if x > threshold else 0, shear))
    volumes = read_cell_volumes(zone, time, case_dir)
    return np.mean(np.array(cells_above_threshold) * np.array(volumes))
