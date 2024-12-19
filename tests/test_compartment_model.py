import unittest

import numpy as np
from of_compartments.cellgrowth.compartment_model import CompartmentModel


class TestCompartmentModel(unittest.TestCase):
    def test_get_flux_array(self):
        concentrations = np.array([1, 2, 1, 0.5])
        F = np.array(
            [[0, 0.1, 0.5, 1.2], [1, 0, 0.1, 0], [0.5, 1, 0, 0.3], [0.3, 0, 1.2, 0]]
        )
        volumes = np.array([1, 2, 1, 1])
        flux = CompartmentModel._get_flux_array(concentrations, F, volumes)
        expected = np.array([0.85, -0.55, -0.5, 0.75])
        np.testing.assert_allclose(flux, expected)
