"""Creates the functions dictionary to put in controlDict for postProcessing"""
import os

import of_compartments.openfoam._function_object_templates as fot
from of_compartments import utils

from . import COMPARTMENT_AVERAGE_FIELDS


def create_function_objects(case_dir, compartment_config, rho, mu, d_ref, p_ref, D):
    compartments = utils.read_compartment_config(compartment_config)

    _write_new_field_control_dict(case_dir, rho, mu, d_ref, p_ref, D)

    _write_compartment_value_control_dict(compartments, case_dir)


def _write_new_field_control_dict(case_dir, rho, mu, d_ref, p_ref, D):
    """Write a set of function objects to create new fields."""
    read_fields_str = f"""
    {fot.read_fields(['V', 'epsilonMean.water', 'alphaMean.air'])}
    {fot.get_tau(rho=rho, mu=mu)}
    {fot.get_kLa(rho=rho, mu=mu, d_ref=d_ref, p_ref=p_ref, D=D)}
    """
    with open(
        os.path.join(case_dir, "system", "controlDict.newFields"), "w"
    ) as output_file:
        output_file.write(fot.HEADER)
        output_file.write(fot.FUNCTIONS_TEMPLATE.format(functions=read_fields_str))


def _write_compartment_value_control_dict(compartments, case_dir):
    """Write function objects to calculate compartment values."""
    functions_str = "#includeFunc writeCellVolumes()\n"
    functions_str += fot.read_fields(["V"] + COMPARTMENT_AVERAGE_FIELDS) + "\n"

    # Calculate volume averaged values and shear stress
    for compartment in compartments:
        functions_str += fot.create_volume_average(
            compartment.id, COMPARTMENT_AVERAGE_FIELDS
        )
        # Write per-cell values for fields to calculate volume fractions
        functions_str += fot.create_tau(compartment.id)
        functions_str += fot.create_volume(compartment.id)

    # Write to new controlDict
    with open(
        os.path.join(case_dir, "system", "controlDict.compartments"), "w"
    ) as output_file:
        output_file.write(fot.HEADER)
        output_file.write(fot.FUNCTIONS_TEMPLATE.format(functions=functions_str))
