# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Unit tests for the signal.py module"""

from datetime import datetime
from pathlib import Path
from time import time
from types import GeneratorType
from unittest.mock import patch, call, PropertyMock, Mock, mock_open

import numpy as np
from freezegun import freeze_time
from pytest import approx, mark, fixture, raises

from spectro_inlets_quantification.config import Config
from spectro_inlets_quantification.exceptions import PeakFitError
from spectro_inlets_quantification.medium import Medium
from spectro_inlets_quantification.peak import Peak
from spectro_inlets_quantification.signal import SignalDict, SignalProcessor

CONFIG = Config()
SIGNAL_MOD = "spectro_inlets_quantification.signal"
FLOAT_RANGE_10 = np.arange(10, dtype=float)


@fixture
def signal_dict():
    return SignalDict(tstamp=47.0, max_depth=2)


@fixture
def signal_processor():
    """Fixture which provides a SignalProcessor"""
    return SignalProcessor(mass_list=["M2", "M44"])


class TestSignalDict:
    """Test the SignalDict class"""

    @freeze_time(datetime(2022, 8, 26, 10, 38, 21, 183143))
    @mark.parametrize("tstamp", [47.0, None])
    @mark.parametrize("max_depth", [47, None])
    @mark.parametrize("verbose", [object(), None])
    def test_init(self, tstamp, max_depth, verbose):
        """Test __init__"""
        init_kwargs = {}
        if tstamp:
            init_kwargs["tstamp"] = tstamp
        if max_depth:
            init_kwargs["max_depth"] = max_depth
        if verbose is not None:
            init_kwargs["verbose"] = verbose

        signal_dict = SignalDict(**init_kwargs)

        assert signal_dict.tstamp == (
            approx(47.0) if tstamp else approx(1661510301.183143)  # Matches frozen time
        )
        assert signal_dict.verbose is (verbose if verbose is not None else False)
        assert signal_dict.max_depth == max_depth
        assert signal_dict.S_vecs == {}
        assert signal_dict.t_vecs == {}

    def test_signals(self, signal_dict):
        """Test the signals method"""
        signal_dict.S_vecs = {
            "M2": np.array([1.23, 4.56, 7.89]),
            "M44": np.array([1.22, 9.88]),
        }
        signals = signal_dict.signals
        assert signals == {"M2": approx(7.89), "M44": approx(9.88)}

    def test_set_signal(self, signal_dict):
        """Test set_signal"""
        assert signal_dict.S_vecs == {}
        assert signal_dict.t_vecs == {}
        signal_dict.set_signal("M2", S_M=6.65, t=101.0)
        assert signal_dict.S_vecs["M2"] == approx(np.array([6.65]))
        assert signal_dict.t_vecs["M2"] == approx(np.array([101.0]))
        signal_dict.set_signal("M2", S_M=7.89, t=105.0)
        assert signal_dict.S_vecs["M2"] == approx(np.array([6.65, 7.89]))
        assert signal_dict.t_vecs["M2"] == approx(np.array([101.0, 105.0]))
        # Add one more to exceed max depth, which will cause the arrays to be the last two
        # calls
        signal_dict.set_signal("M2", S_M=7.89, t=105.0)
        assert signal_dict.S_vecs["M2"] == approx(np.array([7.89, 7.89]))
        assert signal_dict.t_vecs["M2"] == approx(np.array([105.0, 105.0]))

    @mark.parametrize("t", (47.0, None))
    def test_set_signals(self, signal_dict, t):
        """Test set_signals"""
        with patch(SIGNAL_MOD + ".SignalDict.set_signal") as mock_set_signal:
            if t:
                signal_dict.set_signals({"M2": 42, "M44": 47}, t=t)
            else:
                signal_dict.set_signals({"M2": 42, "M44": 47})

            assert mock_set_signal.mock_calls == [
                call("M2", 42, t=t),
                call("M44", 47, t=t),
            ]

    def test_get_signal_non_existing(self, signal_dict):
        """Test get_signal"""
        with raises(KeyError):
            signal_dict.get_signal("non-existing-key")

    @mark.parametrize("tspan", [(2.5, 4.5), None])
    def test_get_signal(self, signal_dict, tspan):
        """Test get_signal"""
        with raises(KeyError):
            signal_dict.get_signal("non-existing-key")

        signal_dict.S_vecs = {"M2": FLOAT_RANGE_10 + 10.0}
        signal_dict.t_vecs = {"M2": FLOAT_RANGE_10}

        t, S_M = signal_dict.get_signal("M2", tspan=tspan)
        if tspan:
            assert t == approx(np.array([3.0, 4.0]))
            assert S_M == approx(np.array([13.0, 14.0]))
        else:
            assert t == approx(FLOAT_RANGE_10)
            assert S_M == approx(FLOAT_RANGE_10 + 10.0)

    @mark.parametrize("tspan", [(2.5, 4.5), None])
    def test_get_signals(self, signal_dict, tspan):
        """Test get_signals"""
        with patch(SIGNAL_MOD + ".SignalDict.get_signal") as mock_get_signal:
            mock_get_signal.side_effect = [3, 4]
            got_signals = signal_dict.get_signals(["M2", "M44"], tspan=tspan)
            assert mock_get_signal.mock_calls == [call("M2", tspan), call("M44", tspan)]
            assert got_signals == {"M2": 3, "M44": 4}

    def test_mass_list(self, signal_dict):
        """Test mass_list"""
        signal_dict.S_vecs = {"M2": None, "M44": None}
        assert signal_dict.mass_list == ["M2", "M44"]

    def test_signals_at(self, signal_dict):
        """Test signals_at"""
        signal_dict.S_vecs = {"M2": FLOAT_RANGE_10, "M44": FLOAT_RANGE_10 + 10.0}
        with raises(IndexError):
            signal_dict.signals_at(50)
        assert signal_dict.signals_at(2) == {"M2": approx(2.0), "M44": approx(12.0)}

    @freeze_time(datetime(2022, 8, 26, 10, 38, 21, 183143))
    @mark.parametrize("N,t", ((3, None), (None, 4.2)))
    def test_get_average_of_last(self, signal_dict, N, t):
        """Test get_average_of_last"""
        t0 = time()
        signal_dict.S_vecs = {
            "M2": FLOAT_RANGE_10 + 10.0,
            "M44": FLOAT_RANGE_10 + 100.0,
        }
        signal_dict.t_vecs = {
            "M2": FLOAT_RANGE_10,
            "M44": FLOAT_RANGE_10,
        }
        signal_dict.tstamp = t0 - 10.0

        S_avg = signal_dict.get_average_of_last(N=N, mass_list=["M2"], t=t)
        if t is not None:
            assert S_avg == {"M2": approx(sum([16.0, 17.0, 18.0, 19.0]) / 4)}
        else:
            assert S_avg == {"M2": approx(sum([17.0, 18.0, 19.0]) / 3)}

    def test_iter(self, signal_dict):
        """Test the __iter__ method"""
        assert isinstance(signal_dict.__iter__(), GeneratorType)
        with patch(
            SIGNAL_MOD + ".SignalDict.mass_list", new_callable=PropertyMock
        ) as mock_mass_list:
            mock_mass_list.return_value = ["mock", "return", "value"]
            collected = []
            for mass in signal_dict:
                collected.append(mass)
            mock_mass_list.assert_called_once()
        assert collected == ["mock", "return", "value"]

    def test_contains(self, signal_dict):
        """Test the __contains__ method"""
        with patch(
            SIGNAL_MOD + ".SignalDict.mass_list", new_callable=PropertyMock
        ) as mock_mass_list:
            mock_mass_list.return_value = ["mock", "return", "value"]
            assert "Not in signal dict" not in signal_dict
            assert "return" in signal_dict

    def test_items(self, signal_dict):
        """Test the items method"""
        assert isinstance(signal_dict.items(), GeneratorType)
        signal_dict.S_vecs = {
            "M2": FLOAT_RANGE_10 + 10.0,
            "M44": FLOAT_RANGE_10 + 100.0,
        }
        with patch(
            SIGNAL_MOD + ".SignalDict.mass_list", new_callable=PropertyMock
        ) as mock_mass_list:
            mock_mass_list.return_value = ["M2", "M44"]
            assert list(signal_dict.items()) == [("M2", 19.0), ("M44", 109.0)]

    @patch(SIGNAL_MOD + ".SignalDict.mass_list", new_callable=PropertyMock)
    @patch(SIGNAL_MOD + ".SignalDict.signals_at")
    @patch(SIGNAL_MOD + ".SignalDict.signals", new_callable=PropertyMock)
    def test_get_item(self, mock_signals, mock_signals_at, mock_mass_list, signal_dict):
        """Test __getitem__"""
        mock_mass_list.return_value = ["M2"]
        mock_signals.return_value = {"M2": 17}
        mock_signals_at.return_value = 19
        signal_dict.t_vecs = {"M2": [1, 2, 3]}
        signal_dict.S_vecs = {"M2": [11, 12, 13]}

        mock_signals.assert_not_called()
        assert signal_dict["M2"] == 17
        mock_signals.assert_called_once()

        mock_signals_at.assert_not_called()
        assert signal_dict[9] == 19
        mock_signals_at.assert_called_once_with(9)

        assert signal_dict["M2-x"] == [1, 2, 3]
        assert signal_dict["M2-y"] == [11, 12, 13]

        with raises(KeyError):
            signal_dict["Not in signal dict"]

    def test_set_item(self, signal_dict):
        """Test __setitem__"""
        with patch(SIGNAL_MOD + ".SignalDict.set_signal") as mock_set_signal:
            signal_dict["some mass"] = 47
            mock_set_signal.assert_called_once_with("some mass", 47)

    @patch(SIGNAL_MOD + ".SignalDict.__contains__")
    def test_clear(self, mock_contains, signal_dict):
        """Test clear"""
        signal_dict.S_vecs = {
            "M2": FLOAT_RANGE_10 + 10.0,
            "M44": FLOAT_RANGE_10 + 100.0,
        }
        signal_dict.t_vecs = {
            "M2": FLOAT_RANGE_10 + 1_000.0,
            "M44": FLOAT_RANGE_10 + 10_000.0,
        }

        mock_contains.return_value = False
        signal_dict.clear("M18")
        assert signal_dict.S_vecs == {
            "M2": approx(FLOAT_RANGE_10 + 10.0),
            "M44": approx(FLOAT_RANGE_10 + 100.0),
        }
        assert signal_dict.t_vecs == {
            "M2": approx(FLOAT_RANGE_10 + 1_000.0),
            "M44": approx(FLOAT_RANGE_10 + 10_000.0),
        }

        mock_contains.return_value = True
        signal_dict.clear("M2")
        assert signal_dict.S_vecs == {"M44": approx(FLOAT_RANGE_10 + 100.0)}
        assert signal_dict.t_vecs == {"M44": approx(FLOAT_RANGE_10 + 10_000.0)}

    @patch(SIGNAL_MOD + ".SignalDict.__iter__")
    def test_clear_all(self, mock_iter, signal_dict):
        """Test clear"""
        signal_dict.S_vecs = {
            "M2": FLOAT_RANGE_10 + 10.0,
            "M44": FLOAT_RANGE_10 + 100.0,
        }
        signal_dict.t_vecs = {
            "M2": FLOAT_RANGE_10 + 1_000.0,
            "M44": FLOAT_RANGE_10 + 10_000.0,
        }

        def iter_mass():
            """Fake iterator over masses"""
            for mass in ["M2", "M44"]:
                yield mass

        mock_iter.return_value = iter_mass()
        signal_dict.clear_all()
        assert signal_dict.S_vecs == {}
        assert signal_dict.t_vecs == {}


class TestSignalProcessor:
    """Unit tests for the SignalProcessor"""

    @mark.parametrize("mass_list", [[1, 2, 3], None])
    @mark.parametrize("peak_type", ["fancy_peak", None])
    @mark.parametrize("max_depth", [123, None])
    @mark.parametrize("nonlin_coeff", [(1.23, 4.56), None])
    @mark.parametrize("signal_dict", [Mock(), None])
    @mark.parametrize("tstamp", [47.0, None])
    @mark.parametrize("verbose", [True, None])
    def test_init(
        self,
        mass_list,
        peak_type,
        max_depth,
        nonlin_coeff,
        signal_dict,
        tstamp,
        verbose,
    ):
        """Test the __init__"""
        init_kwargs = {}
        for argname, value in locals().items():
            if argname in ("self", "init_kwargs"):
                continue
            if value is not None:
                init_kwargs[argname] = value

        with patch(SIGNAL_MOD + ".SignalDict") as mock_signal_dict_class:
            returned_mock_signal_dict = Mock()
            mock_signal_dict_class.return_value = returned_mock_signal_dict
            signal_processors = SignalProcessor(**init_kwargs)

        expected_mass_list = [1, 2, 3] if mass_list else []
        assert signal_processors.mass_list == expected_mass_list

        expected_peak_type = peak_type if peak_type else "gauss"
        assert signal_processors.peak_type == expected_peak_type

        expected_nonlin_coeff = nonlin_coeff if nonlin_coeff else (0, 0)
        assert signal_processors.nonlin_coeff == expected_nonlin_coeff

        # A signal dict is created in __init__ if not provided
        expected_verbose = verbose if verbose is not None else False
        if signal_dict:
            assert signal_processors.signal_dict is signal_dict
        else:
            assert signal_processors.signal_dict is returned_mock_signal_dict
            expected_max_depth = max_depth if max_depth else 1000
            mock_signal_dict_class.assert_called_once_with(
                max_depth=expected_max_depth, tstamp=tstamp, verbose=expected_verbose
            )

        assert signal_processors.verbose is expected_verbose
        assert isinstance(signal_processors.medium, Medium)

    @patch("builtins.open", new_callable=mock_open, read_data='{"key": "value"}')
    @mark.parametrize("file_name", ["some filename", None])
    @mark.parametrize("proc_dir", [Path("home"), None])
    def test_load(self, mock_open, file_name, proc_dir, reset_singletons):
        """Test the load class method"""
        init_kwargs = {}
        if file_name:
            init_kwargs["file_name"] = file_name
        if proc_dir:
            init_kwargs["proc_dir"] = proc_dir

        # expected_proc_dir = proc_dir if proc_dir else CONFIG.processor_directories[0]
        expected_filename = Path((file_name if file_name else "Sorens final processor") + ".yml")

        with patch("spectro_inlets_quantification.config.Config.get_best_data_file") as get_best:
            get_best.return_value = Path("Best")
            with patch(SIGNAL_MOD + ".SignalProcessor.__init__") as signal_processor_init:
                signal_processor_init.return_value = None
                SignalProcessor.load(**init_kwargs)
        get_best.assert_called_once_with(
            data_file_type="processors", filepath=expected_filename, override_source_dir=proc_dir
        )

        mock_open.assert_called_once_with(Path("Best"), "r")
        signal_processor_init.assert_called_once_with(key="value")

    def test_signals(self, signal_processor):
        """Test signals property"""
        signal_processor.signal_dict = 47
        assert signal_processor.signals == 47

    def test_pvac(self, signal_processor):
        """Test p_vac"""
        signal_processor.medium.p_vac = 56
        assert signal_processor.p_vac == 56

    def test_tstamp(self, signal_processor):
        """Test tstamp"""
        signal_processor.signal_dict.tstamp = 89
        assert signal_processor.tstamp == 89

    def test_peak_class(self, signal_processor):
        """Test the PeakClass property"""
        with patch(SIGNAL_MOD + ".PEAK_CLASSES", {"peaktype": 11}):
            signal_processor.peak_type = "peaktype"
            assert signal_processor.PeakClass == 11

            signal_processor.peak_type = "non existent peak type"
            with raises(KeyError):
                signal_processor.PeakClass

        signal_processor.peak_type = None
        assert signal_processor.PeakClass is Peak

    @mark.parametrize("p_vac", [1.23e-4, None])
    def test_cals_nonlinear_factor(self, p_vac, signal_processor):
        """Test calc_nonlinear_factor"""
        signal_processor.nonlin_coeff = (47.47, 1.23)
        p_vac = p_vac if p_vac else signal_processor.p_vac
        p_hat = p_vac / 1e-4
        non_linear_factor = 1 + 47.47 * p_hat + 1.23 * p_hat**2
        assert non_linear_factor == approx(signal_processor.calc_nonlinear_factor(p_vac=p_vac))

    @mark.parametrize("p_vac", [47.47, None])
    def test_correct_y(self, p_vac, signal_processor):
        """Test correct_y"""
        # NOTE To expand if background correction is implemented
        with patch.object(signal_processor, "calc_nonlinear_factor") as mock_calc:
            mock_calc.return_value = 47.0
            assert signal_processor.correct_y(
                FLOAT_RANGE_10,
                FLOAT_RANGE_10 + 123.4,
                p_vac=p_vac,
            ) == approx((FLOAT_RANGE_10 + 123.4) / 47.0)
            mock_calc.assert_called_once_with(p_vac=p_vac)

    @freeze_time(datetime(2022, 8, 26, 10, 38, 21, 183143))
    @patch(SIGNAL_MOD + ".Peak")
    @mark.parametrize("t", [47.0, None])
    @mark.parametrize("Mspan", ((10.8, 11.2), None))
    @mark.parametrize("fit_successful", [True, False])
    def test_make_peak(self, mock_peak, t, Mspan, fit_successful, signal_processor):
        """Test make_peak"""
        x = np.linspace(start=10.0, stop=12.0, num=21)
        y = x * 123.4
        if Mspan:
            x_masked = np.linspace(Mspan[0], Mspan[1], 5)
        else:
            x_masked = np.linspace(10.5, 11.5, 11)
        y_masked = x_masked * 123.4

        # Mock peaks
        mock_gauss_peak = Mock()
        PEAK_CLASSES = {"gauss": mock_gauss_peak}

        # Create a fake tstamp
        signal_processor.signal_dict.tstamp = 1661510301.183143 - 100.0
        if t is None:
            t_expected = 100.0
        else:
            t_expected = t

        fake_peak_result = object()
        if fit_successful:
            mock_gauss_peak.return_value = fake_peak_result
        else:
            mock_gauss_peak.side_effect = PeakFitError
            mock_peak.return_value = fake_peak_result

        with patch.object(signal_processor, "correct_y") as mock_correct:
            with patch(SIGNAL_MOD + ".PEAK_CLASSES", PEAK_CLASSES):
                mock_correct.return_value = y_masked
                peak = signal_processor.make_peak(x, y, mass="M11", Mspan=Mspan, t=t)
            mock_correct.assert_called_once_with(approx(x_masked), approx(y_masked), p_vac=None)

        if fit_successful:
            mock_gauss_peak.assert_called_once_with(
                x=approx(x_masked), y=approx(y_masked), t=t_expected
            )
        else:
            mock_peak.assert_called_once_with(
                x=approx(x_masked), y=approx(y_masked), t=t_expected, error=True
            )

        assert peak is fake_peak_result

    @mark.parametrize("Mspan_in", ((10.8, 11.2), None))
    @mark.parametrize("fit_width_in", (1.2, None))
    def test_calc_signal(self, Mspan_in, fit_width_in, signal_processor):
        """Test calc_signal"""
        x = np.linspace(10, 12, 20)
        y = x + 123.4

        def fake_make_peak(x, y, mass, Mspan=None, fit_width=1.0):
            """Fake make peak"""
            assert Mspan == Mspan_in
            assert fit_width == fit_width
            peak = Peak(x, y)
            signal_processor.mass_peaks["M11"] = peak
            return peak

        with patch.object(signal_processor, "make_peak", fake_make_peak):
            S_M = signal_processor.calc_signal(
                x=x, y=y, mass="M11", Mspan=Mspan_in, fit_width=fit_width_in
            )

        # Since we just return a normal peak, where the signal is max(y)
        assert S_M == approx(12.0 + 123.4)

    @mark.parametrize("fit_width", (1.2, None))
    def test_calc_signals(self, fit_width, signal_processor):
        """Test calc_signals"""
        x = FLOAT_RANGE_10
        y = FLOAT_RANGE_10 + 123.4
        signal_processor.signal_dict = 67
        with patch.object(signal_processor, "calc_signal") as mock_calc_signal:
            kwargs = {"fit_width": fit_width} if fit_width else {}
            return_value = signal_processor.calc_signals(x, y, ["M2", "M44"], **kwargs)
            expected_fitwidth = fit_width if fit_width else 1.0
            assert mock_calc_signal.mock_calls == [
                call(x=x, y=y, mass="M2", fit_width=expected_fitwidth),
                call(x=x, y=y, mass="M44", fit_width=expected_fitwidth),
            ]
        assert return_value == 67

    @mark.parametrize("mass_list", (["M13", "M5"], None))
    def test_get_average_of_last(self, mass_list, signal_processor):
        """Test get_average_of_last"""
        expected_mass_list = mass_list if mass_list else signal_processor.mass_list
        with patch.object(signal_processor.signal_dict, "get_average_of_last") as mock_average:
            mock_average.return_value = 47
            return_value = signal_processor.get_average_of_last(123, mass_list=mass_list)
        assert return_value == 47
        mock_average.assert_called_once_with(123, expected_mass_list)

    def test_getattr_(self, signal_processor):
        """Test __getattr__"""
        # Should raise on empty mass_peaks
        signal_processor.mass_peaks = {}
        with raises(AttributeError):
            signal_processor.prop

        class Fake:
            def __init__(self, prop):
                self.prop = prop

        signal_processor.mass_peaks = {"M2": Fake(42), "M44": Fake(47)}
        # Should map the proper from each peak
        assert signal_processor.prop == {"M2": 42, "M44": 47}

        # And still raise for an unknown property
        with raises(AttributeError):
            signal_processor.prop2
