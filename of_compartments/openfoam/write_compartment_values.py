"""Writes compartment values from openFOAM postprocessing to .csv files"""
import os

import of_compartments.openfoam.compartment_data_reader as cdr
import pandas as pd
from of_compartments import utils


def write_compartment_values(
    case_dir, compartment_config, output_dir, time, high_shear_threshold
):
    """Write all compartment average and interface values to .csv"""
    compartment_df, interface_df = get_compartment_values(
        time, case_dir, compartment_config, high_shear_threshold
    )

    os.system(f"mkdir -p {output_dir}")
    compartment_output_path = os.path.join(output_dir, "compartment_values.csv")
    interface_output_path = os.path.join(output_dir, "interface_values.csv")
    compartment_df.to_csv(compartment_output_path, index=False)
    interface_df.to_csv(interface_output_path, index=False)


def get_compartment_values(time, case_dir, compartment_config, high_shear_threshold):
    """Returns dictionary for each field mapping compartment ids to values"""
    compartments = utils.read_compartment_config(compartment_config)
    compartment_ids = utils.get_ids(compartments)
    top_compartments = list(utils.get_top_compartments_by_id(compartment_ids))

    # Compartment averages
    compartment_dict = {
        "compartment": compartment_ids,
        "gas_holdup": [
            cdr.read_gas_holdup(id, time, case_dir) for id in compartment_ids
        ],
        "kLa": [cdr.read_kLa(id, time, case_dir) for id in compartment_ids],
        "epsilon": [cdr.read_epsilon(id, time, case_dir) for id in compartment_ids],
        "volume": [cdr.read_volume(id, time, case_dir) for id in compartment_ids],
        "top_gas_flux": cdr.read_top_gas_flux(compartment_ids, case_dir),
        "tau": [cdr.read_tau_average(id, time, case_dir) for id in compartment_ids],
        "high_tau_fraction": [
            cdr.read_tau_threshold_fraction(id, time, high_shear_threshold, case_dir)
            for id in compartment_ids
        ],
        "is_top": [id in top_compartments for id in compartment_ids],
    }
    compartment_df = pd.DataFrame.from_dict(compartment_dict)

    # Interface values
    flux_matrix = cdr.read_flow(compartment_ids, case_dir, "water")
    print("\t- Correcting net flow")
    corrected_flow = utils.correct_flow(flux_matrix)

    compartment_idx_tuples = [
        (i, j) for i in range(len(compartment_ids)) for j in range(len(compartment_ids))
    ]
    interface_dict = {
        "compartment_src": [compartment_ids[c[0]] for c in compartment_idx_tuples],
        "compartment_dest": [compartment_ids[c[1]] for c in compartment_idx_tuples],
        "corrected_flow": [corrected_flow[c[0], c[1]] for c in compartment_idx_tuples],
    }

    interface_df = pd.DataFrame.from_dict(interface_dict)

    return compartment_df, interface_df
