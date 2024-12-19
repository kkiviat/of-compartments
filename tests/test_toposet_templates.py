import os
import unittest

import of_compartments.openfoam._toposet_templates as tst
import of_compartments.openfoam.create_compartment_toposet as ts
from of_compartments.utils import Compartment


class TestToposetTemplates(unittest.TestCase):
    def test_write_toposet_dict(self):
        case_dir = os.path.join("tests", "data", "test_case")
        ts.create_compartment_toposet(
            case_dir, os.path.join(case_dir, "compartment_config")
        )

        with open(os.path.join(case_dir, "system", "topoSetDict.expected"), "r") as f:
            expected = f.readlines()
        with open(os.path.join(case_dir, "system", "topoSetDict"), "r") as f:
            result = f.readlines()

        self.assertEqual(result, expected)

    def test_get_compartment_string(self):
        compartment = Compartment(
            bottom_height=1,
            top_height=1.5,
            inner_radius=1.1,
            outer_radius=1.2,
            id="h1r2",
        )
        result = tst.get_compartment_string(compartment)
        expected = """
{
    name    h1r2;
    type    cellSet;
    action  new;

    source  cylinderAnnulusToCell;
    p1 (0.0 1 0.0);
    p2 (0.0 1.5 0.0);
    innerRadius 1.1;
    outerRadius 1.2;
}

{
    name        h1r2_zone;
    type        cellZoneSet;
    action      new;

    source      setToCellZone;
    sourceInfo
    {
        set   h1r2;
    }
}
"""
        self.assertEqual(result, expected)

    def test_get_boundary_face_string(self):
        """
        Test that it creates the topoSetDict string for boundary faces
        """
        compartment1 = Compartment(
            bottom_height=1,
            top_height=1.5,
            inner_radius=1.1,
            outer_radius=1.2,
            id="h1r2",
        )
        compartment2 = Compartment(
            bottom_height=0,
            top_height=1.5,
            inner_radius=0,
            outer_radius=1.1,
            id="h0r0",
        )
        result = tst.get_boundary_face_string(compartment1, compartment2)
        self.maxDiff = None
        expected = """
{
    name        boundary_h1r2_h0r0;
    type        faceSet;
    action      new;

    source      cellToFace;
    sourceInfo {
      set h1r2;
      option all;
    }
}
{
    name        boundary_h1r2_h0r0;
    type        faceSet;
    action      subset;

    source      cellToFace;
    sourceInfo {
      set h0r0;
      option all;
    }
}
// Convert faceSet to faceZone
{
    name        boundary_h1r2_h0r0_zone;
    type        faceZoneSet;
    action      new;

    source setsToFaceZone;
    sourceInfo
    {
        faceSet boundary_h1r2_h0r0;
        cellSet h1r2;
    }
}
"""
        self.assertEqual(result, expected)

    def test_get_top_face_string(self):
        """
        Test that it gets the topoSetDict entry for the top face zone
        """
        compartment = Compartment(
            bottom_height=1,
            top_height=1.5,
            inner_radius=1.1,
            outer_radius=1.2,
            id="h1r2",
        )
        result = tst.get_top_face_string(compartment)
        expected = """
{
    name        boundary_h1r2_top;
    type        faceSet;
    action      new;

    source      cellToFace;
    sourceInfo {
      set h1r2;
      option all;
    }
}
{
    name    boundary_h1r2_top;
    type    faceSet;
    action  subset;

    source  cylinderToFace;
    p1 (0.0 1.4985 0.0);
    p2 (0.0 3.0 0.0);
    radius 1.2;
}
{
    name        boundary_h1r2_top;
    type        faceSet;
    action      subset;

    source normalToFace;
    sourceInfo
    {
        normal (0 1 0);     // Vector
        cos     0.1;       // Tolerance (max cos of angle)
    }
}
// Convert faceSet to faceZone
{
    name        boundary_h1r2_top_zone;
    type        faceZoneSet;
    action      new;

    source setsToFaceZone;
    sourceInfo
    {
        faceSet boundary_h1r2_top;
        cellSet h1r2;
    }
}
"""
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
