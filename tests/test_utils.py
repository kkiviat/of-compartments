import unittest

import numpy as np
from of_compartments import utils
from of_compartments.utils import Compartment


class TestUtils(unittest.TestCase):
    def test_get_height_index_from_compartment_id(self):
        self.assertEqual(utils.get_height_index_from_compartment_id("h1r2"), 1)
        self.assertEqual(utils.get_height_index_from_compartment_id("h10r20"), 10)

    def test_get_radius_index_from_compartment_id(self):
        self.assertEqual(utils.get_radius_index_from_compartment_id("h1r2"), 2)
        self.assertEqual(utils.get_radius_index_from_compartment_id("h10r20"), 20)

    def test_get_ids(self):
        compartments = []
        compartments.append(self.make_compartment_id("h0r0"))
        compartments.append(self.make_compartment_id("h0r1"))
        compartments.append(self.make_compartment_id("h1r0"))
        compartments.append(self.make_compartment_id("h1r1"))
        self.assertEqual(utils.get_ids(compartments), ["h0r0", "h0r1", "h1r0", "h1r1"])

    def test_get_compartment_lower_upper_pairs(self):
        compartments = []
        compartments.append(self.make_compartment_id("h0r0"))
        compartments.append(self.make_compartment_id("h0r1"))
        compartments.append(self.make_compartment_id("h0r2"))
        compartments.append(self.make_compartment_id("h1r0"))
        compartments.append(self.make_compartment_id("h1r1"))
        compartments.append(self.make_compartment_id("h1r2"))
        compartments.append(self.make_compartment_id("h2r0"))
        compartments.append(self.make_compartment_id("h2r1"))
        compartments.append(self.make_compartment_id("h2r2"))
        self.assertCountEqual(
            utils.get_compartment_lower_upper_pairs(compartments),
            [
                ("h0r0", "h1r0"),
                ("h0r1", "h1r1"),
                ("h0r2", "h1r2"),
                ("h1r0", "h2r0"),
                ("h1r1", "h2r1"),
                ("h1r2", "h2r2"),
            ],
        )

    def test_order_compartment(self):
        self.assertEqual(
            utils.order_compartments("h0r1", "h1r0"),
            ("h0r1", "h1r0"),
        )
        self.assertEqual(
            utils.order_compartments("h1r0", "h0r1"),
            ("h0r1", "h1r0"),
        )
        self.assertEqual(
            utils.order_compartments("h1r0", "h1r1"),
            ("h1r0", "h1r1"),
        )
        self.assertEqual(
            utils.order_compartments("h1r1", "h1r0"),
            ("h1r0", "h1r1"),
        )

    def test_correct_flow(self):
        flow = np.array([[0, 1.01], [1.02, 0]])
        corrected_flow = utils.correct_flow(flow)
        print(corrected_flow)
        self.assertAlmostEqual(corrected_flow[0, 1], corrected_flow[1, 0], places=4)
        self.assertAlmostEqual(flow[0, 1], corrected_flow[0, 1], places=1)
        self.assertAlmostEqual(flow[1, 0], corrected_flow[1, 0], places=1)

    def make_compartment_id(self, id, is_top=False):
        return Compartment(
            bottom_height=0,
            top_height=1,
            inner_radius=0,
            outer_radius=1,
            id=id,
        )


if __name__ == "__main__":
    unittest.main()
