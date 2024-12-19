"""Clone the given case to the target location.

Includes the latest (parallel) time step, reconstructs
latest time, and removes processor directories
"""
import argparse
import os
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        "-c",
        help="source case directory",
        required=True,
    )
    parser.add_argument(
        "--target",
        "-t",
        help="name of directory for copied case",
        required=True,
    )

    args = parser.parse_args()
    case_dir = args.case
    target = args.target
    subprocess.check_call(
        f"pyFoamCloneCase.py --no-pyfoam --latest-time --parallel {case_dir} {target}",
        shell=True,
    )
    subprocess.check_call(f"reconstructPar -latestTime -case {target}", shell=True)
    subprocess.check_call(f"rm -r {target}/processor*", shell=True)
    os.mknod(os.path.join(target, ".cloned"))

    # Overwrite controlDict with a generic one in case there is something weird
    # or location-specific about the original
    with open(os.path.join(target, "system", "controlDict"), "w") as f:
        f.write(
            """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  10
     \\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format          ascii;
    class           dictionary;
    location        "system";
    object          controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     multiphaseEulerFoam;
startFrom       latestTime;
startTime       0;
stopAt          endTime;
endTime         90;
deltaT          1;
writeControl    adjustableRunTime;
writeInterval   1;
purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   7;
runTimeModifiable yes;
            """
        )


if __name__ == "__main__":
    main()
