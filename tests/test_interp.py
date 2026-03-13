"""Tests for pymeteo.interp module."""

import numpy as np
import pytest
import pymeteo.interp as interp


class TestLinear:
    """Test linear interpolation."""

    def test_exact_point(self):
        dim = np.array([0.0, 1.0, 2.0, 3.0])
        var = np.array([10.0, 20.0, 30.0, 40.0])
        result = interp.linear(dim, var, 1.0)
        assert abs(result - 20.0) < 1e-10

    def test_midpoint(self):
        dim = np.array([0.0, 1.0, 2.0, 3.0])
        var = np.array([10.0, 20.0, 30.0, 40.0])
        result = interp.linear(dim, var, 1.5)
        assert abs(result - 25.0) < 1e-10

    def test_quarter_point(self):
        dim = np.array([0.0, 4.0])
        var = np.array([0.0, 100.0])
        result = interp.linear(dim, var, 1.0)
        assert abs(result - 25.0) < 1e-10

    def test_handles_masked_values(self):
        dim = np.array([0.0, 1.0, np.nan, 3.0])
        var = np.array([0.0, 10.0, 20.0, 30.0])
        # Should still work, skipping NaN in dim
        result = interp.linear(dim, var, 2.0)
        assert np.isfinite(result)

    def test_dimension_mismatch_raises(self):
        dim = np.array([0.0, 1.0, 2.0])
        var = np.array([10.0, 20.0])
        with pytest.raises(Exception, match="Dimensions"):
            interp.linear(dim, var, 1.0)


class TestInterpHeight:
    """Test height interpolation to a pressure level."""

    def test_exact_pressure_level(self):
        z = np.array([0.0, 1000.0, 2000.0, 3000.0])
        p = np.array([100000.0, 90000.0, 80000.0, 70000.0])
        result = interp.interp_height(z, p, 90000.0)
        assert abs(result - 1000.0) < 1e-6

    def test_below_surface_returns_zero(self):
        z = np.array([0.0, 1000.0, 2000.0])
        p = np.array([100000.0, 90000.0, 80000.0])
        result = interp.interp_height(z, p, 105000.0)
        assert result == 0

    def test_above_top_returns_neg1(self):
        z = np.array([0.0, 1000.0, 2000.0])
        p = np.array([100000.0, 90000.0, 80000.0])
        result = interp.interp_height(z, p, 50000.0)
        assert result == -1

    def test_interpolated_value_between_levels(self):
        z = np.array([0.0, 1000.0, 2000.0, 3000.0])
        p = np.array([100000.0, 90000.0, 80000.0, 70000.0])
        result = interp.interp_height(z, p, 85000.0)
        # Should be between 1000 and 2000
        assert 1000.0 < result < 2000.0


class TestInterpPressure:
    """Test pressure interpolation to a height level."""

    def test_exact_height_level(self):
        p = np.array([100000.0, 90000.0, 80000.0, 70000.0])
        z = np.array([0.0, 1000.0, 2000.0, 3000.0])
        result = interp.interp_pressure(p, z, 1000.0)
        assert abs(result - 90000.0) < 1e-3

    def test_interpolated_value(self):
        p = np.array([100000.0, 90000.0, 80000.0, 70000.0])
        z = np.array([0.0, 1000.0, 2000.0, 3000.0])
        result = interp.interp_pressure(p, z, 500.0)
        # Should be between 100000 and 90000
        assert 90000.0 < result < 100000.0

    def test_above_top_returns_none(self):
        p = np.array([100000.0, 90000.0, 80000.0])
        z = np.array([0.0, 1000.0, 2000.0])
        result = interp.interp_pressure(p, z, 5000.0)
        assert result is None
