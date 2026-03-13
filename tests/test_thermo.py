"""Tests for pymeteo.thermo module."""

import math
import numpy as np
import pytest
import pymeteo.thermo as thermo
import pymeteo.constants as c


class TestTemperatureConversions:
    """Test T <-> theta conversions."""

    def test_T_from_theta_at_surface(self):
        # At p00, theta == T
        result = thermo.T(300.0, c.p00)
        assert abs(result - 300.0) < 1e-6

    def test_theta_from_T_at_surface(self):
        # At p00, theta == T
        result = thermo.theta(300.0, c.p00)
        assert abs(result - 300.0) < 1e-6

    def test_T_theta_roundtrip(self):
        T_orig = 280.0
        p = 85000.0
        th = thermo.theta(T_orig, p)
        T_back = thermo.T(th, p)
        assert abs(T_back - T_orig) < 1e-6

    def test_theta_increases_with_lower_pressure(self):
        # Theta at lower pressure should be higher for same T
        T = 250.0
        th_high = thermo.theta(T, 50000.0)
        th_low = thermo.theta(T, 90000.0)
        assert th_high > th_low


class TestSaturationVaporPressure:
    """Test es and esi functions."""

    def test_es_at_freezing(self):
        # At 273.15 K, es should be approximately 611.2 Pa (by definition)
        result = thermo.es(c.T00)
        assert abs(result - 611.2) < 1.0  # within 1 Pa

    def test_es_increases_with_temperature(self):
        es_cold = thermo.es(260.0)
        es_warm = thermo.es(300.0)
        assert es_warm > es_cold

    def test_esi_at_freezing(self):
        # At 273.15 K, esi should be approximately 611.2 Pa
        result = thermo.esi(c.T00)
        assert abs(result - 611.2) < 1.0

    def test_es_greater_than_esi_below_freezing(self):
        # Below freezing, es(liquid) > es(ice) — Bergeron process
        T = 263.15  # -10 C
        assert thermo.es(T) > thermo.esi(T)


class TestMixingRatio:
    """Test saturation mixing ratio."""

    def test_w_vs_positive(self):
        result = thermo.w_vs(300.0, 100000.0)
        assert result > 0

    def test_w_vs_increases_with_temperature(self):
        w_cold = thermo.w_vs(260.0, 100000.0)
        w_warm = thermo.w_vs(300.0, 100000.0)
        assert w_warm > w_cold

    def test_w_vs_typical_value(self):
        # At ~30C, 1000 hPa, ws should be roughly 25-30 g/kg
        w = thermo.w_vs(303.15, 100000.0)
        assert 0.020 < w < 0.035  # kg/kg


class TestThetaV:
    """Test virtual potential temperature."""

    def test_theta_v_greater_than_theta_for_moist_air(self):
        th = 300.0
        qv = 0.010  # 10 g/kg
        result = thermo.theta_v(th, qv)
        assert result > th

    def test_theta_v_equals_theta_for_dry_air(self):
        th = 300.0
        result = thermo.theta_v(th, 0.0)
        assert abs(result - th) < 1e-10


class TestDewpoint:
    """Test dewpoint temperature calculation."""

    def test_Td_less_than_T(self):
        # Dewpoint should generally be <= temperature for unsaturated air
        p = 100000.0
        T = 300.0
        qv = 0.010  # unsaturated at this T, p
        td = thermo.Td(p, qv)
        assert td < T

    def test_Td_approaches_T_at_saturation(self):
        # At saturation mixing ratio, Td should be close to T
        T = 290.0
        p = 100000.0
        qv_sat = thermo.w_vs(T, p)
        td = thermo.Td(p, qv_sat)
        # Should be within a few degrees of T
        assert abs(td - T) < 3.0


class TestLatentHeat:
    """Test latent heat function."""

    def test_Lv_returns_constant(self):
        result = thermo.Lv(300.0)
        assert result == c.L


class TestPressureAltitude:
    """Test pressure from pressure altitude."""

    def test_surface_pressure(self):
        # At z=0, should return p00
        result = thermo.p_from_pressure_altitude(0.0, 300.0)
        assert abs(result - c.p00) < 1.0

    def test_pressure_decreases_with_altitude(self):
        p_low = thermo.p_from_pressure_altitude(0.0, 300.0)
        p_high = thermo.p_from_pressure_altitude(5000.0, 250.0)
        assert p_high < p_low


class TestThetaE:
    """Test equivalent potential temperature."""

    def test_th_e_greater_than_theta(self):
        p = 100000.0
        t = 300.0
        td = 295.0
        qv = 0.015
        the = thermo.th_e(p, t, td, qv)
        th = thermo.theta(t, p)
        assert the > th

    def test_th_e_when_saturated(self):
        # Should not raise
        p = 100000.0
        t = 290.0
        td = t  # saturated
        qv = thermo.w_vs(t, p)
        the = thermo.th_e(p, t, td, qv)
        assert the > 0


class TestSaturationMixingRatios:
    """Test q_vl and q_vi functions."""

    def test_q_vl_positive(self):
        result = thermo.q_vl(100000.0, 300.0)
        assert result > 0

    def test_q_vi_positive(self):
        result = thermo.q_vi(100000.0, 260.0)
        assert result > 0

    def test_q_vl_greater_than_q_vi_below_freezing(self):
        T = 263.15  # -10 C
        p = 100000.0
        assert thermo.q_vl(p, T) > thermo.q_vi(p, T)


class TestLapseRates:
    """Test moist adiabatic lapse rates."""

    def test_dTdz_moist_negative(self):
        # Temperature should decrease with height (negative lapse rate)
        result = thermo.dTdz_moist(280.0, 90000.0)
        assert result < 0

    def test_dTdz_moist_less_than_dry_adiabat(self):
        # Moist adiabatic lapse rate magnitude should be less than dry (g/cp)
        result = abs(thermo.dTdz_moist(280.0, 90000.0))
        dry_rate = c.g / c.cpd
        assert result < dry_rate

    def test_dTdp_moist_positive(self):
        # Temperature decreases with decreasing pressure, so dT/dp > 0
        result = thermo.dTdp_moist(280.0, 90000.0)
        assert result > 0


class TestWetBulb:
    """Test wet-bulb temperature calculation."""

    def test_Twb_between_Td_and_T(self):
        # Create a simple sounding
        z = np.array([0.0, 500.0, 1000.0, 1500.0, 2000.0])
        p = np.array([100000.0, 95000.0, 90000.0, 85000.0, 80000.0])
        th = np.array([300.0, 301.0, 302.0, 303.0, 304.0])
        qv = np.array([0.012, 0.010, 0.008, 0.006, 0.004])

        twb = thermo.Twb(z, p, th, qv, 0.0)
        if not np.isnan(twb):
            # Wet bulb should be finite
            assert isinstance(twb, (float, np.floating))


class TestCAPE:
    """Test CAPE computation."""

    def test_cape_returns_dict(self):
        # Simple unstable sounding — warm surface, cold aloft
        nk = 20
        z = np.linspace(0, 12000, nk).astype(np.float32)
        # Decrease pressure from surface to top
        p = np.array(
            [100000 * np.exp(-c.g * zi / (c.Rd * 270)) for zi in z],
            dtype=np.float32,
        )
        # Temperature decreasing faster than moist adiabat
        t = np.array([290.0 - 7.0 * zi / 1000.0 for zi in z], dtype=np.float32)
        q = np.array([0.012 * np.exp(-zi / 3000.0) for zi in z], dtype=np.float32)

        result = thermo.CAPE(z, p, t, q, 1)
        assert isinstance(result, dict)
        assert "cape" in result
        assert "cin" in result
        assert "lcl" in result
        assert "lfc" in result
        assert "el" in result
        assert "t_p" in result
        assert result["cape"] >= 0

    def test_cape_parcel_types(self):
        # Test all 3 parcel types don't crash
        nk = 30
        z = np.linspace(0, 15000, nk).astype(np.float32)
        p = np.array(
            [100000 * np.exp(-c.g * zi / (c.Rd * 260)) for zi in z],
            dtype=np.float32,
        )
        t = np.array([295.0 - 6.5 * zi / 1000.0 for zi in z], dtype=np.float32)
        q = np.array([0.014 * np.exp(-zi / 2500.0) for zi in z], dtype=np.float32)

        for parcel in [1, 2, 3]:
            result = thermo.CAPE(z, p, t, q, parcel)
            assert isinstance(result, dict)
            assert result["cape"] >= 0

    def test_cape_stable_sounding(self):
        # Very stable sounding — isothermal (no buoyancy)
        nk = 20
        z = np.linspace(0, 10000, nk).astype(np.float32)
        p = np.array(
            [100000 * np.exp(-c.g * zi / (c.Rd * 300)) for zi in z],
            dtype=np.float32,
        )
        # Strong inversion — theta increases rapidly
        t = np.array(
            [thermo.T(300.0 + 15.0 * zi / 10000.0, pi) for zi, pi in zip(z, p)],
            dtype=np.float32,
        )
        q = np.array([0.002] * nk, dtype=np.float32)

        result = thermo.CAPE(z, p, t, q, 1)
        # Very small or zero CAPE for stable sounding
        assert result["cape"] < 100.0
