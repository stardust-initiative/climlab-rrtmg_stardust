"""Synthetic algebra checks, runnable with the Python standard unittest runner."""
import unittest

import numpy as np

from offline_scattering import Slab, demonstration, join, reflected, solve_interfaces, stack


class SlabChecks(unittest.TestCase):
    def test_dense_solve_independent_of_composition(self):
        rng = np.random.default_rng(42)
        for n in (1, 2, 4, 8):
            layers = [Slab(*(np.array([.03, .12, .7, .1, .65, .08, .5])
                              * rng.uniform(.8, 1., 7))) for _ in range(n)]
            combined = stack(layers)
            beam, down, up = solve_interfaces(layers, .9, .2, .15)
            self.assertAlmostEqual(up[0], combined.reflect_beam*.9
                                   + combined.reflect_down*.2 + combined.transmit_up*.15, places=14)
            self.assertAlmostEqual(down[-1], combined.source_beam*.9
                                   + combined.transmit_down*.2 + combined.reflect_up*.15, places=14)
            self.assertAlmostEqual(beam[-1], combined.transmit_beam*.9, places=14)

    def test_identity_and_associativity(self):
        a = Slab(.03, .12, .7, .1, .65, .08, .5)
        b = Slab(.04, .09, .75, .12, .7, .06, .6)
        self.assertEqual(join(Slab.identity(), a), a)
        self.assertEqual(join(a, Slab.identity()), a)
        x, y = join(join(a, b), a), join(a, join(b, a))
        np.testing.assert_allclose(list(vars(x).values()), list(vars(y).values()), atol=1e-14, rtol=0)

    def test_lambertian_cancellation_and_two_response_sensitivity(self):
        result = demonstration()
        self.assertLess(result['interface_solve_max_error'], 1e-14)
        self.assertLess(result['lambertian_formula_error'], 1e-14)
        self.assertLess(abs(result['lambertian_beam_split_difference']), 1e-14)
        self.assertGreater(abs(result['two_response_beam_split_difference']), .01)

    def test_transparent_absorbing_and_conservative_limits(self):
        self.assertAlmostEqual(reflected(Slab.identity(), .3, .4), .3)
        absorbing = Slab.symmetric_cell(0, 0, .5, .5, .5)
        self.assertAlmostEqual(reflected(absorbing, .3, .4), .075)
        conservative = Slab.symmetric_cell(.2, .2, .8, .8, .5)
        b, d, u = solve_interfaces([conservative, conservative])
        self.assertTrue(np.all(b >= 0) and np.all(d >= 0) and np.all(u >= 0))
        self.assertAlmostEqual(u[0]+d[-1]+b[-1], 1., places=14)

    def test_invalid_inputs_not_clipped(self):
        with self.assertRaises(ValueError):
            Slab(float('nan'), 0, 1, 0, 1, 0, 1)
        with self.assertRaises(ValueError):
            reflected(Slab.identity(), 1.1, .3)
        with self.assertRaises(ValueError):
            reflected(Slab(0, 0, 1, 0, 1, 1, 1), 1, 1)


if __name__ == '__main__':
    unittest.main()
