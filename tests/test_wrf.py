"""Tests for pymeteo.wrf module."""

import numpy as np
import pytest
import pymeteo.wrf as wrf


class TestLlToIj:
    """Test latitude/longitude to grid index conversion."""

    def test_reference_point_returns_origin(self):
        # At the reference point, (i, j) should be near (0, 0) or (1, 1)
        # depending on convention. For Lambert:
        result = wrf.ll_to_ij(
            map_proj=1,
            truelat1=30.0,
            truelat2=60.0,
            stand_lon=-98.0,
            dx=12000.0,
            dy=12000.0,
            ref_lat=40.0,
            ref_lon=-98.0,
            lat=40.0,
            lon=-98.0,
        )
        assert result is not None
        i, j = result
        # Should be very close to origin (0 or 1)
        assert abs(i) <= 1
        assert abs(j) <= 1

    def test_lambert_conformal(self):
        result = wrf.ll_to_ij(
            map_proj=1,
            truelat1=30.0,
            truelat2=60.0,
            stand_lon=-98.0,
            dx=12000.0,
            dy=12000.0,
            ref_lat=40.0,
            ref_lon=-98.0,
            lat=41.0,
            lon=-97.0,
        )
        assert result is not None
        assert len(result) == 2

    def test_polar_stereographic(self):
        result = wrf.ll_to_ij(
            map_proj=2,
            truelat1=60.0,
            truelat2=90.0,
            stand_lon=-80.0,
            dx=32000.0,
            dy=32000.0,
            ref_lat=60.0,
            ref_lon=-80.0,
            lat=62.0,
            lon=-78.0,
        )
        assert result is not None
        assert len(result) == 2

    def test_mercator(self):
        result = wrf.ll_to_ij(
            map_proj=3,
            truelat1=0.0,
            truelat2=0.0,
            stand_lon=0.0,
            dx=30000.0,
            dy=30000.0,
            ref_lat=0.0,
            ref_lon=0.0,
            lat=1.0,
            lon=1.0,
        )
        assert result is not None
        assert len(result) == 2

    def test_unsupported_projection_returns_none(self):
        # map_proj=99 is unsupported
        result = wrf.ll_to_ij(
            map_proj=99,
            truelat1=30.0,
            truelat2=60.0,
            stand_lon=-98.0,
            dx=12000.0,
            dy=12000.0,
            ref_lat=40.0,
            ref_lon=-98.0,
            lat=40.0,
            lon=-98.0,
        )
        assert result is None

    def test_offset_from_reference(self):
        # Moving east should increase i
        ref_result = wrf.ll_to_ij(
            map_proj=1,
            truelat1=30.0,
            truelat2=60.0,
            stand_lon=-98.0,
            dx=12000.0,
            dy=12000.0,
            ref_lat=40.0,
            ref_lon=-98.0,
            lat=40.0,
            lon=-98.0,
        )
        east_result = wrf.ll_to_ij(
            map_proj=1,
            truelat1=30.0,
            truelat2=60.0,
            stand_lon=-98.0,
            dx=12000.0,
            dy=12000.0,
            ref_lat=40.0,
            ref_lon=-98.0,
            lat=40.0,
            lon=-95.0,  # 3 degrees east
        )
        assert east_result[0] > ref_result[0]
