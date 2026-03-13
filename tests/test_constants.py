"""Tests for pymeteo.constants module."""

import pymeteo.constants as c


class TestConstants:
    """Verify that physical constants have expected values and relationships."""

    def test_gas_constants(self):
        assert c.Rd == 287.04
        assert c.Rv == 461.5
        assert abs(c.epsilon - c.Rd / c.Rv) < 1e-10

    def test_gravity(self):
        assert c.gravity == 9.81
        assert c.g == 9.81

    def test_latent_heat(self):
        assert c.L == 2.501e6
        assert c.xlv == c.L

    def test_heat_capacities(self):
        assert c.cp == 1005.7
        assert c.cpd == 1005.7
        assert c.cpv == 1875.0
        assert c.cv == 718.0

    def test_reference_values(self):
        assert c.p00 == 100000.0
        assert c.T00 == 273.15

    def test_derived_kappa(self):
        expected_kappa = (c.cp - c.cv) / c.cp
        assert abs(c.kappa - expected_kappa) < 1e-10

    def test_derived_kappa_d(self):
        expected = c.Rd / c.cp
        assert abs(c.kappa_d - expected) < 1e-10

    def test_derived_rp00(self):
        assert abs(c.rp00 - 1.0 / c.p00) < 1e-20

    def test_derived_reps(self):
        assert abs(c.reps - c.Rv / c.Rd) < 1e-10

    def test_conversion_factors(self):
        assert c.m2km == 0.001
        assert c.km2m == 1000.0

    def test_converge(self):
        assert c.converge == 0.0002
