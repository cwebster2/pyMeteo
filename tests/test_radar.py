"""Tests for pymeteo.radar module."""

import numpy as np
import pytest
import pymeteo.radar as radar


class TestRadarLevels:
    """Test radar level functions."""

    def test_cmap_radar_levels_full_range(self):
        levels = radar.cmap_radar_levels_full()
        assert levels[0] == pytest.approx(-25.0)
        assert levels[-1] == pytest.approx(74.5)
        # Should have (75 - (-25)) / 0.5 = 200 entries
        assert len(levels) == 200

    def test_cmap_radar_levels_range(self):
        levels = radar.cmap_radar_levels()
        assert levels[0] == pytest.approx(0.0)
        assert levels[-1] == pytest.approx(74.5)
        # Should have 75 / 0.5 = 150 entries
        assert len(levels) == 150


class TestRadarColors:
    """Test radar color palette."""

    def test_cmap_radar_colors_count(self):
        colors = radar.cmap_radar_colors()
        assert len(colors) == 15

    def test_cmap_radar_colors_are_hex(self):
        colors = radar.cmap_radar_colors()
        for color in colors:
            assert color.startswith("#")
            assert len(color) == 7


class TestContourLevels:
    """Test contour level function."""

    def test_cmap_contour_levels(self):
        levels = radar.cmap_contour_levels()
        assert levels == [30, 40, 50, 60, 65]


class TestColormapRegistration:
    """Test that colormaps are registered."""

    def test_radar_cmap_registered(self):
        import matplotlib

        try:
            cmap = matplotlib.colormaps["pymeteo_radar"]
            assert cmap is not None
        except (KeyError, AttributeError):
            # Older matplotlib
            import matplotlib.cm as cm

            cmap = cm.get_cmap("pymeteo_radar")
            assert cmap is not None

    def test_radar_full_cmap_registered(self):
        import matplotlib

        try:
            cmap = matplotlib.colormaps["pymeteo_radar_full"]
            assert cmap is not None
        except (KeyError, AttributeError):
            import matplotlib.cm as cm

            cmap = cm.get_cmap("pymeteo_radar_full")
            assert cmap is not None
