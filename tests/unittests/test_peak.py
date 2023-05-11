# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Unit tests for the peak.py module"""

from unittest.mock import patch, PropertyMock

import numpy as np
import pytest
from pytest import mark, fixture, approx

from spectro_inlets_quantification.constants import STANDARD_WIDTH
from spectro_inlets_quantification.exceptions import PeakFitError
from spectro_inlets_quantification.peak import Peak, gauss_fun, GaussPeak


@fixture
def peak():
    """Fixture that provides an initialized Peak object"""
    # Simulate a peak with half a sine
    x = np.linspace(
        0.0, np.pi, 101
    )  # Odd number of points to make sure there is an x for max y
    return Peak(x=x, y=np.sin(x) + 1.0, t=47.0, bg=1.0)


@fixture
def gauss_peak(peak):
    """A fixture that provides a gauss peak"""
    kwargs = {
        "center": 10.0,
        "sigma": 2.3,
        "x": peak.x,
        "y": peak.y,
        "t": 10.0,
        "tolerance": 10000.0,
    }
    return GaussPeak(**kwargs)


PEAK_MOD = "spectro_inlets_quantification.peak"


class TestPeak:
    """Test the peak class"""

    @mark.parametrize("bg", [None, 42.0])
    def test_init(self, bg):
        """Test the Peak.__init__"""
        x = np.arange(10)
        y = np.arange(10, 20)
        t = 47.0
        error = object()
        peak = Peak(x=x, y=y, t=t, bg=bg, error=error)
        assert peak.x == approx(x)
        assert peak.y_raw == approx(y)
        assert peak.t == approx(t)
        if bg:
            assert peak.y == approx(y - bg)
            assert peak.bg == approx(bg)
        else:
            assert peak.y == approx(y)
            assert peak.bg is None

        for cache_name in ("height", "center", "width", "y_fit", "x_fit", "integral"):
            assert getattr(peak, f"_{cache_name}") is None

        assert peak.error is error

    @mark.parametrize("background", (1.23, [1.23] * 101, np.array([1.23] * 101)))
    @mark.parametrize(
        "x", (np.linspace(0.0, np.pi, 51), list(np.linspace(0.0, np.pi, 51)))
    )
    def test_background_function(self, peak, background, x):
        """Test background function"""
        peak.bg = background
        calculated_background = peak.background_function(x)
        assert calculated_background == approx(np.array([1.23] * 51))

    @mark.parametrize("property_", ("height", "center", "width", "integral"))
    def test_calculated_cached_properties(self, property_, peak):
        """Test the calculated, cached properties; height, center, width, integral"""
        cache_property_name = f"_{property_}"
        calc_method_name = f"calc_{property_}"
        with patch.object(peak, calc_method_name) as calc_mock:
            calc_mock.return_value = 34.5

            # Get value once
            value = getattr(peak, property_)
            assert value == approx(calc_mock.return_value)
            calc_mock.assert_called_once()
            assert getattr(peak, cache_property_name) == approx(calc_mock.return_value)
            # Get value once more
            value2 = getattr(peak, property_)
            assert value == approx(value2)
            # The calc method should still only have been called once due to caching
            calc_mock.assert_called_once()

    def test_y_fit(self, peak):
        """Test y_fit"""
        peak._y_fit = object()
        assert peak.y_fit == peak.y_fit

    def test_x_fit(self, peak):
        """Test x_fit"""
        x_fit = peak.x_fit
        assert x_fit == approx(np.arange(0.0, np.pi, 0.025))
        second_x_fit = peak.x_fit
        # Test caching
        assert x_fit is second_x_fit

    def test_calc_height(self, peak):
        """Test calc_height, (simple max of y)"""
        assert peak.calc_height() == approx(1.0)

    def test_calc_center(self, peak):
        """Test calc_center (simple x of max y)"""
        # sine has max at pi / 2
        assert peak.calc_center() == approx(np.pi / 2)

    def test_calc_width(self, peak):
        """Test calc_width"""
        assert peak.calc_width() == approx(2.8902652)

    def test_calc_integral(self, peak):
        """Test calc_integral"""
        assert peak.calc_integral() == approx(2.0, abs=0.01)

    @mark.parametrize("error", (None, object()))
    def test_signal(self, error, peak):
        """Test signal"""
        peak.error = error
        peak._height = 7.0
        # 0 if peak has error, else height
        if error:
            assert peak.signal == approx(0.0)
        else:
            assert peak.signal == approx(7.0)


def test_gauss_fun():
    """Test the gauss function"""
    for height in np.arange(1e-12, 1e-11, 1e-12):
        for sigma in np.arange(0.1, 10.1, 1):
            for center in np.arange(-5, 205, 20):
                x = np.arange(-5, 205, 10)
                assert gauss_fun(
                    x, center=center, sigma=sigma, height=height
                ) == approx(height * np.exp(-((x - center) ** 2) / (2 * sigma**2)))


class TestGaussPeak:
    """Test the GaussPeak class"""

    @patch(PEAK_MOD + ".Peak.__init__")
    @patch(PEAK_MOD + ".GaussPeak.fit_gauss")
    @mark.parametrize("tolerance", (30000.3, None))
    def test_init(self, mock_fit, mock_init, peak, tolerance):
        """Test the __init__ method"""
        kwargs = {
            "center": 10.0,
            "sigma": 2.3,
            "x": peak.x,
            "y": peak.y,
            "t": 10.0,
        }
        if tolerance:
            kwargs["tolerance"] = tolerance
        gauss_peak = GaussPeak(**kwargs)
        if tolerance:
            assert gauss_peak.tolerance == approx(tolerance)
        else:
            assert gauss_peak.tolerance == approx(0.2)
        mock_fit.assert_called_once_with(center=10.0, sigma=2.3)
        mock_init.assert_called_once_with(peak.x, peak.y, t=10.0)
        assert gauss_peak._fwhm is None
        assert gauss_peak.sigma is None

    @patch(
        PEAK_MOD + ".GaussPeak.relative_square_error",
        new_callable=PropertyMock,
    )
    @patch(PEAK_MOD + ".GaussPeak.background_function")
    @patch(PEAK_MOD + ".GaussPeak.x_fit", new_callable=PropertyMock)
    @patch(PEAK_MOD + ".curve_fit")
    @mark.parametrize("center", [10.0, None])
    @mark.parametrize("sigma", [2.3, None])
    @mark.parametrize("relative_square_error", [0.1, 10000.0])
    def test_fit_gauss(
        self,
        mock_curve_fit,
        mock_x_fit,
        mock_background_function,
        mock_relative_square_error,
        peak,
        center,
        sigma,
        relative_square_error,
    ):
        """Test fit gauss peak"""
        # This test function is long and difficult to reason about, however I couldn't at the time figure out a good way
        # to break it up, so I suggest reading it in parallel with the method under test \Kenneth

        mock_x_fit.return_value = np.arange(10)
        mock_background_function.return_value = np.array([0.0] * 10)
        mock_relative_square_error.return_value = relative_square_error

        number_of_argument_provided = len(
            [arg for arg in (center, sigma) if arg is not None]
        )
        # First argument returned by curve_fit are the fitted parameters, filled in below
        if number_of_argument_provided == 2:
            # Both center and sigma are provided, so we are only fitting the height, which is
            # then the only parameter returned
            mock_curve_fit.return_value = ([23], 45)
        elif number_of_argument_provided == 1:
            mock_curve_fit.return_value = ([12, 33], 45)
        else:  # 0 parameters given
            mock_curve_fit.return_value = ([56, 67, 78], 45)

        # Instantiate and make sure bad fits raises
        kwargs = {
            "center": center,
            "sigma": sigma,
            "x": peak.x,
            "y": peak.y,
            "t": 10.0,
        }
        if relative_square_error > 0.2:
            with pytest.raises(PeakFitError):
                GaussPeak(**kwargs)
            return
        else:
            gauss_peak = GaussPeak(**kwargs)

        guess_c = (peak.x[-1] + peak.x[0]) / 2
        guess_s = STANDARD_WIDTH / 3
        guess_h = max(peak.y)

        # Extract arguments to curve_fit and test those
        mock_curve_fit.assert_called_once()
        call_args = mock_curve_fit.call_args
        fit_function, x, y = call_args[0]
        p0 = call_args[1]["p0"]

        assert x == approx(peak.x)
        assert y == approx(peak.y)

        with patch(PEAK_MOD + ".gauss_fun") as mock_gauss_fun:
            if number_of_argument_provided == 2:
                assert fit_function.__name__ == "gauss_i_height_only"
                assert p0 == guess_h
                height = 23
                # Test calling though the reduced test function
                fit_function(6, 8)
                mock_gauss_fun.assert_called_once_with(
                    6, center=center, sigma=sigma, height=8
                )
            elif number_of_argument_provided == 1:
                if center is not None:
                    assert fit_function.__name__ == "gauss_i_sigma_height"
                    assert p0 == [guess_s, guess_h]
                    sigma, height = 12, 33
                    # Test calling though the reduced test function
                    fit_function(6, 8, 9)
                    mock_gauss_fun.assert_called_once_with(
                        6, center=center, sigma=8, height=9
                    )
                else:
                    assert fit_function.__name__ == "gauss_i_center_height"
                    assert p0 == [guess_c, guess_h]
                    center, height = 12, 33
                    # Test calling though the reduced test function
                    fit_function(6, 8, 9)
                    mock_gauss_fun.assert_called_once_with(
                        6, center=8, sigma=sigma, height=9
                    )
            else:
                assert fit_function.__name__ == "gauss_i_all"
                assert p0 == [guess_c, guess_s, guess_h]
                center, sigma, height = 56, 67, 78
                # Test calling though the reduced test function
                fit_function(6, 8, 9, 10)
                mock_gauss_fun.assert_called_once_with(6, center=8, sigma=9, height=10)

        # Make sure the provided test results of the fit were written to the appropriate properties
        assert gauss_peak.height == approx(height)
        assert gauss_peak.center == approx(center)
        assert gauss_peak.sigma == approx(sigma)

        assert gauss_peak.y_fit == approx(
            gauss_fun(np.arange(10), center, sigma, height)
        )
        assert gauss_peak.integral == approx(np.sqrt(2 * np.pi) * height * sigma)

    def test_gauss_fit_curve_fit_raise(self, peak):
        """Test the case where the curve_fit in gauss_fit raises"""
        kwargs = {
            "x": peak.x,
            "y": peak.y,
            "t": 10.0,
        }
        with patch(PEAK_MOD + ".curve_fit") as mock_curve_fit:
            mock_curve_fit.side_effect = RuntimeError
            with pytest.raises(PeakFitError):
                GaussPeak(**kwargs)

    @patch(PEAK_MOD + ".gauss_fun")
    def test_y_of_x(self, mock_gauss_fun, gauss_peak):
        """Test y_of_x"""
        x = np.arange(10)
        gauss_peak.y_of_x(x)
        mock_gauss_fun.assert_called_once()
        args, kwargs = mock_gauss_fun.call_args
        assert x == approx(args[0])
        # Look at gauss_peak fixture for values
        assert kwargs["center"] == approx(10.0)
        assert kwargs["sigma"] == approx(2.3)
        # TODO the test value below is merely read of from current implementation, not predicted
        assert kwargs["height"] == approx(
            74.771983
        )  # Rather poor fit of gauss peak to sine

    @patch(PEAK_MOD + ".GaussPeak.y_of_x")
    @patch(PEAK_MOD + ".GaussPeak.background_function")
    def test_y_raw_of_x(self, mock_background_function, mock_y_of_x, gauss_peak):
        """Test y_raw_of_x"""
        range_10 = np.arange(10)
        mock_y_of_x.return_value = range_10 + 11.0
        mock_background_function.return_value = range_10 + 13.0
        result = gauss_peak.y_raw_of_x(range_10)
        assert mock_y_of_x.call_args[0][0] == approx(range_10)
        assert mock_background_function.call_args[0][0] == approx(range_10)
        assert result == approx(range_10 * 2 + 24.0)

    @patch(PEAK_MOD + ".GaussPeak.y_of_x")
    def test_relative_square_error(self, mock_y_of_x, gauss_peak):
        """Test relative_square_error"""
        range_10 = np.arange(10)
        gauss_peak.x = range_10
        gauss_peak.y = np.sin(range_10) + 0.1
        mock_y_of_x.return_value = np.sin(range_10) + 0.3
        result = gauss_peak.relative_square_error
        assert mock_y_of_x.call_args[0][0] == approx(range_10)
        error = np.array([0.3 - 0.1] * 10)  # The sines cancel out
        assert result == approx(error.dot(error) / gauss_peak.y.dot(gauss_peak.y))

    def test_fwhm(self, gauss_peak):
        """Test fwhm (full width at half maximum) property"""
        gauss_peak.sigma = 34.5
        result = gauss_peak.fwhm
        assert result == approx(2 * np.sqrt(2 * np.log(2)) * 34.5)

    def test_width(self, gauss_peak):
        """Test the width property"""
        gauss_peak.sigma = 34.5
        assert gauss_peak.width == approx(34.5)
