"""Create a topoSetDict to define compartments"""
import os
import subprocess

import of_compartments.openfoam._toposet_templates as tst
from of_compartments import utils


def create_compartment_toposet(case_dir, compartment_config):
    compartments = utils.read_compartment_config(compartment_config)

    # clear old sets and zones
    subprocess.run(
        f"rm {os.path.join(case_dir, 'constant', 'polyMesh', 'sets', '*')}",
        shell=True,
        check=False,
    )
    subprocess.run(
        f"rm {os.path.join(case_dir, 'constant', 'polyMesh', 'cellZones')}",
        shell=True,
        check=False,
    )
    subprocess.run(
        f"rm {os.path.join(case_dir, 'constant', 'polyMesh', 'faceZones')}",
        shell=True,
        check=False,
    )

    # Remove MRFProperties in case it exists since we're deleting any MRF zones
    MRF_file = os.path.join(case_dir, "constant", "MRFProperties")
    if os.path.isfile(MRF_file):
        subprocess.run(f"rm {MRF_file}", shell=True, check=False)

    # Write topoSetDict
    with open(os.path.join(case_dir, "system", "topoSetDict.compartments"), "w") as f:
        f.write(tst.HEADER)

        # Define compartments
        for compartment in compartments:
            f.write(tst.get_compartment_string(compartment))

        # Define boundary faces
        for i in range(len(compartments)):
            for j in range(i + 1, len(compartments)):
                f.write(tst.get_boundary_face_string(compartments[i], compartments[j]))

        # top zones
        for compartment in utils.get_top_compartments(compartments):
            f.write(tst.get_top_face_string(compartment))

        f.write(tst.FOOTER)
