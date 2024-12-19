import os
import unittest

import of_compartments.openfoam.create_function_objects as cfo


class TestCompartmentFunctionObjects(unittest.TestCase):
    def test_write_control_dict(self):
        case_dir = os.path.join("tests", "data", "test_case")
        cfo.create_function_objects(
            case_dir,
            os.path.join(case_dir, "compartment_config"),
            rho=993,
            mu=0.0007,
            d_ref=0.004,
            p_ref=1e5,
            D=2e-9,
        )

        with open(
            os.path.join(case_dir, "system", "controlDict.newFields.expected"), "r"
        ) as f:
            expected = f.readlines()
        with open(os.path.join(case_dir, "system", "controlDict.newFields"), "r") as f:
            result = f.readlines()

        self.assertEqual(result, expected)

        with open(
            os.path.join(case_dir, "system", "controlDict.compartments.expected"), "r"
        ) as f:
            expected = f.readlines()
        with open(
            os.path.join(case_dir, "system", "controlDict.compartments"), "r"
        ) as f:
            result = f.readlines()

        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
