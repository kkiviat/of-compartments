import unittest

import numpy as np
import of_compartments.cellgrowth.damage_models as dm
from of_compartments import utils
from of_compartments.cellgrowth.compartment_model import CompartmentData


class TestDamageModels(unittest.TestCase):
    def test_eddyDamageModel_lakhotia_papoutsakis(self):
        compartment_data = self.make_compartment_data(1, 2, epsilon=np.array([1, 2]))
        eddyDamageModel = dm.eddyDamageModel_lakhotia_papoutsakis(
            B=1, k_c=10, E_0=0.5, nu=1e-6
        )
        q1, q2 = eddyDamageModel(compartment_data)
        np.testing.assert_allclose(q1, 254.2789028)
        np.testing.assert_allclose(
            q2,
            403.76738260,
        )

    def test_shearDamageModel_constant(self):
        compartment_data = self.make_compartment_data(
            1, 2, high_shear_fraction=np.array([0.1, 0.8])
        )
        shearDamageModel = dm.shearDamageModel_constant(0.2)
        q1, q2 = shearDamageModel(compartment_data)
        np.testing.assert_allclose(q1, 0.02)
        np.testing.assert_allclose(q2, 0.16)

    def test_bubbleDamageModel_wall2017(self):
        compartment_ids = ["h0r0", "h0r1", "h1r0", "h1r1"]
        volumes = [0.1, 1, 1, 0.5]
        compartment_data = self.make_compartment_data(
            2,
            2,
            volume=volumes,
            compartment_ids=compartment_ids,
            top_gas_flux=np.array([0.0007, 1, 0.01, 0.004]),
        )
        bubbleDamageModel = dm.bubbleDamageModel_walls2017(0.002, "low")
        q1, q2, q3, q4 = bubbleDamageModel(compartment_data)

        # Bulk compartments should have 0 death rate
        self.assertEqual(q1, 0)
        self.assertEqual(q2, 0)
        # Top compartments should have positive death rate
        self.assertGreater(q3, 0)
        self.assertGreater(q4, 0)

    def test_bubbleDamageModel_cherry1992(self):
        compartment_ids = ["h0r0", "h0r1", "h1r0", "h1r1"]
        volumes = np.array([0.1, 1, 1, 0.5])
        compartment_data = self.make_compartment_data(
            2,
            2,
            volume=volumes,
            compartment_ids=compartment_ids,
            top_gas_flux=np.array([0.0007, 1, 0.01, 0.004]),
        )
        bubbleDamageModel = dm.bubbleDamageModel_cherry1992(0.002)
        q1, q2, q3, q4 = bubbleDamageModel(compartment_data)

        # Bulk compartments should have 0 death rate
        self.assertEqual(q1, 0)
        self.assertEqual(q2, 0)
        # Top compartments should have positive death rate
        self.assertGreater(q3, 0)
        self.assertGreater(q4, 0)

    def make_compartment_data(
        self,
        n_radial,
        n_vertical,
        compartment_ids=None,
        is_top=None,
        gas_holdup=None,
        high_shear_fraction=None,
        kLa=None,
        volume=None,
        top_gas_flux=None,
        epsilon=None,
    ):
        if compartment_ids is None:
            compartment_ids = [
                utils.get_zone_id(h, r)
                for h in range(n_vertical)
                for r in range(n_radial)
            ]
        num_compartments = n_radial * n_vertical
        if is_top is None:
            top_compartments = utils.get_top_compartments_by_id(compartment_ids)
            is_top = np.array([id in top_compartments for id in compartment_ids])
            print(is_top)
        if gas_holdup is None:
            gas_holdup = np.zeros(num_compartments)
        if high_shear_fraction is None:
            high_shear_fraction = np.zeros(num_compartments)
        if kLa is None:
            kLa = np.zeros(num_compartments)
        if volume is None:
            volume = np.zeros(num_compartments)
        if top_gas_flux is None:
            top_gas_flux = np.zeros(num_compartments)
        if epsilon is None:
            epsilon = np.zeros(num_compartments)
        compartment_data = CompartmentData(
            ids=compartment_ids,
            gas_holdup=gas_holdup,
            high_shear_fraction=high_shear_fraction,
            kLa=kLa,
            volumes=volume,
            top_gas_flux=top_gas_flux,
            epsilon=epsilon,
            is_top=is_top,
        )
        return compartment_data


if __name__ == "__main__":
    unittest.main()
