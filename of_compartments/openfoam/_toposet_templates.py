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
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      topoSetDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

actions
(
"""

FOOTER = r"""
);

// ************************************************************************* //
"""


FACE_INTERSECTION = """
{{
    name        boundary_{zone1}_{zone2};
    type        faceSet;
    action      new;

    source      cellToFace;
    sourceInfo {{
      set {zone1};
      option all;
    }}
}}
{{
    name        boundary_{zone1}_{zone2};
    type        faceSet;
    action      subset;

    source      cellToFace;
    sourceInfo {{
      set {zone2};
      option all;
    }}
}}
// Convert faceSet to faceZone
{{
    name        boundary_{zone1}_{zone2}_zone;
    type        faceZoneSet;
    action      new;

    source setsToFaceZone;
    sourceInfo
    {{
        faceSet boundary_{zone1}_{zone2};
        cellSet {zone1};
    }}
}}
"""

TOP_BOUNDARY_FACE_ZONE = """
{{
    name        boundary_{zone}_top;
    type        faceSet;
    action      new;

    source      cellToFace;
    sourceInfo {{
      set {zone};
      option all;
    }}
}}
{{
    name    boundary_{zone}_top;
    type    faceSet;
    action  subset;

    source  cylinderToFace;
    p1 (0.0 {height_bottom} 0.0);
    p2 (0.0 {height_top} 0.0);
    radius {radius};
}}
{{
    name        boundary_{zone}_top;
    type        faceSet;
    action      subset;

    source normalToFace;
    sourceInfo
    {{
        normal (0 1 0);     // Vector
        cos     0.1;       // Tolerance (max cos of angle)
    }}
}}
// Convert faceSet to faceZone
{{
    name        boundary_{zone}_top_zone;
    type        faceZoneSet;
    action      new;

    source setsToFaceZone;
    sourceInfo
    {{
        faceSet boundary_{zone}_top;
        cellSet {zone};
    }}
}}
"""

CELL_SET = """
{{
    name    {zone};
    type    cellSet;
    action  new;

    source  cylinderAnnulusToCell;
    p1 (0.0 {height_bottom} 0.0);
    p2 (0.0 {height_top} 0.0);
    innerRadius {inner_radius};
    outerRadius {outer_radius};
}}
"""

CELL_SET_TO_ZONE = """
{{
    name        {zone}_zone;
    type        cellZoneSet;
    action      new;

    source      setToCellZone;
    sourceInfo
    {{
        set   {zone};
    }}
}}
"""


def get_compartment_string(compartment):
    """Returns string to put in topoSetDict to define cell set / zone for compartment"""
    return CELL_SET.format(
        zone=compartment.id,
        height_bottom=compartment.bottom_height,
        height_top=compartment.top_height,
        inner_radius=compartment.inner_radius,
        outer_radius=compartment.outer_radius,
    ) + CELL_SET_TO_ZONE.format(zone=compartment.id)


def get_boundary_face_string(compartment1, compartment2):
    return FACE_INTERSECTION.format(zone1=compartment1.id, zone2=compartment2.id)


def get_top_face_string(compartment):
    return TOP_BOUNDARY_FACE_ZONE.format(
        zone=compartment.id,
        height_bottom=compartment.top_height * 0.999,
        height_top=compartment.top_height * 2,
        radius=compartment.outer_radius,
    )
