import os
import subprocess

from of_compartments import utils

HEADER = r"""
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  10
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
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


def read_fields(fields):
    return f"""
readFields
{{
    functionObjectLibs ( "libfieldFunctionObjects.so" );
    type readFields;
    fields ({" ".join(fields)});
}}
        """


def get_tau(rho, mu):
    return f"""
tau_s
{{
    functionObjectLibs ( "libutilityFunctionObjects.so" );
    type            coded;

    codeOptions
    #{{
        -I$(LIB_SRC)/meshTools/lnInclude
    #}};

    codeExecute
    #{{
        const volScalarField& eps
            (
            mesh().lookupObject<volScalarField>("epsilonMean.water")
            );

        volScalarField tau_s(2.5 * {mu} * pow(eps / (6 * ({mu}/{rho})), 0.5));
        tau_s.rename("tau_s");
        tau_s.write();
    #}};
}}
    """


def get_kLa(rho, mu, d_ref, p_ref, D):
    return f"""
kLa
{{
    functionObjectLibs ( "libutilityFunctionObjects.so" );
    type            coded;

    codeOptions
    #{{
        -I$(LIB_SRC)/meshTools/lnInclude
    #}};

    codeExecute
    #{{
        const volScalarField& eps
            (
            mesh().lookupObject<volScalarField>("epsilonMean.water")
            );
        const volScalarField& alphaAir
            (
            mesh().lookupObject<volScalarField>("alphaMean.air")
            );
        const volScalarField& p
            (
            mesh().lookupObject<volScalarField>("p")
            );

        volScalarField b_d({d_ref} * pow({p_ref} / p, 1/3));
        volScalarField a(6 * alphaAir / {d_ref});
        volScalarField kL(3600 * 2 * pow({D}/M_PI, 0.5) * pow(eps * {rho}/{mu}, 0.25));
        volScalarField kLa(kL * a);
        kLa.rename("kLa");
        kLa.write();
    #}};
}}
    """


FUNCTIONS_TEMPLATE = """
functions
{{
{functions}
}}
"""

_BOUNDARY_FORMAT = "boundary_{boundary_name}_zone"


def get_boundary_empty(zone1, zone2, case_dir):
    "Check if there are actually any faces in the intersection of these face sets"
    boundary_name = _BOUNDARY_FORMAT.format(
        boundary_name=utils.create_boundary_name(zone1, zone2)
    )
    file_name = os.path.join(case_dir, "constant", "polyMesh", "sets", boundary_name)

    if not os.path.isfile(file_name):
        return True

    # Read 18th line in sets file to get number of faces
    num_faces = int(subprocess.check_output(f"sed '18q;d' {file_name}", shell=True))

    return num_faces == 0


def create_volume_average(zone, average_fields):
    # Calculate fields averaged over the volume
    _VOL_AVG_TEMPLATE = """
    volAvg_{zone}
    {{
        type            volFieldValue;
        libs            ("libfieldFunctionObjects.so");

        writeFields     false;
        writeToFile     true;
        writeControl    writeTime;

        regionType      cellZone;
        name            {zone}_zone;
        operation       volAverage;

        fields          ({fields});
    }}
    """
    return _VOL_AVG_TEMPLATE.format(zone=zone, fields=" ".join(average_fields))


def _create_vol_field_value(name, zone, fields):
    return f"""
    {name}_{zone}
    {{
        type            volFieldValue;
        libs            ("libfieldFunctionObjects.so");

        writeFields     true;
        writeControl    writeTime;

        regionType      cellZone;
        name            {zone}_zone;
        operation       none;

        fields          ({" ".join(fields)});
    }}
    """


def create_tau(zone):
    return _create_vol_field_value("tau", zone, ["tau_s"])


def create_interface_area(boundary_name):
    """Compute a surface field arbitrarily just to get the area"""
    return f"""
    area_{boundary_name}
    {{
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");

        writeFields     true;
        writeControl    writeTime;

        surfaceFormat   none;
        regionType      faceZone;
        name            boundary_{boundary_name}_zone;
        operation       sum;

        fields          (phi.water);
    }}
    """


def create_volume(zone):
    _VOLUME_TEMPLATE = """
    volumes_{zone}
    {{
        type            volFieldValue;
        libs            ("libfieldFunctionObjects.so");

        writeFields     true;
        writeToFile     false;
        writeControl    writeTime;

        regionType      cellZone;
        name            {zone}_zone;
        operation       none;

        fields          (V);
    }}
    """
    return _VOLUME_TEMPLATE.format(zone=zone)
