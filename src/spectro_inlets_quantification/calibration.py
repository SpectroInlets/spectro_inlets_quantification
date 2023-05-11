# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Master module for handling results of calibration experiment

Classes here inherit from classes in sensitivity.

Equivalently initiating a `SensitivityList` from a ``sf_list`` of SensitivityFactors,
the way to make a calibration during dev-time is to add CalPoints to a `cal_list` and
then use that to initiate a `Calibration`.
"""
from pathlib import Path
from typing import Dict, List, Union, Optional, Any, cast, TYPE_CHECKING

import attr
import json
import time
from .constants import STANDARD_MOL_COLORS
from .custom_types import JSONVALUE, PATHLIKE, MOLLIST, MASSLIST, MOL, MASS
from .tools import TODAY, mass_to_setting, mass_to_M, make_axis
from .molecule import MoleculeDict, Molecule
from .sensitivity import (  # noqa
    SensitivityFactor,
    SensitivityUnion,
    SensitivityList,
    SensitivityMatrix,
    SensitivityFit,
    SENSITIVITYLIST_FILTER_TYPE,
)
from .config import Config
from .chip import Chip

if TYPE_CHECKING:
    from matplotlib import pyplot


CONFIG = Config()


@attr.s
class CalPoint(SensitivityFactor):
    """Represents one result of a calibration experiment

    Attributes inherited from SensitivityFactor:

    Attributes:
        mol (str): The name of molecule
        cal_type (str): The type of calibration, i.e. "internal", "semi", "external", etc
        F (dict): ``{mass: F_i_M}`` where ``mass`` is a str and ``F_i_M`` is the sensitivity
            factor in C/mol for mol at mass

    Calibration-specific attributes:

    Attributes:
        precision (float): The relative rms error during the calibration measurement
        background_signal (float): The average background signal at mass during the
            calibration experiment
        background_std (float): The standard deviation of the background signal at mass
            during the calibration experiment. This is used to calculate the
            detection limit (see `calc_detection_limit` and `Calibration.prints_report`)
        description(string): Description of the calibration experiment
        date(string): Date of the calibration experiment
        setup(string): Setup on which the calibration experiment was done
        internal_conditions (dict): Conditions in the vacuum chamber and mass
            spectrometer during calibration, e.g. E_ion, V_CEM, J_fill, roughing pump...
        external_conditions (dict): Conditions in and above the chip during calibration,
            e.g. T, p, carrier gas, solvent, etc
    """

    precision: Optional[float] = attr.ib(default=None)
    background_signal: Optional[float] = attr.ib(default=None)
    background_std: Optional[float] = attr.ib(default=None)
    description: Optional[str] = attr.ib(default=None)
    date: str = attr.ib(default=TODAY)
    setup: Optional[str] = attr.ib(default=None)
    internal_conditions: Optional[Dict[str, JSONVALUE]] = attr.ib(default=None)
    external_conditions: Optional[Dict[str, JSONVALUE]] = attr.ib(default=None)

    def as_dict(self) -> Dict[str, Union[float, str, Dict[str, JSONVALUE]]]:
        """Return the dictionary representation of the CalPoint"""
        self_as_dict = cast(
            Dict[str, Union[float, str, Dict[str, JSONVALUE]]],
            super().as_dict(),
        )
        self_as_dict.update(
            precision=self.precision,
            background_signal=self.background_signal,
            background_std=self.background_std,
            description=self.description,
            date=self.date,
            setup=self.setup,
            internal_conditions=self.internal_conditions,
            external_conditions=self.external_conditions,
        )
        return self_as_dict

    def copy(self) -> "CalPoint":
        """Return a new CalPoint cloning self"""
        self_as_dict = self.as_dict()
        for key, value in self_as_dict.items():
            if isinstance(value, dict):
                self_as_dict[key] = value.copy()
        return CalPoint(**self_as_dict)  # type: ignore

    def calc_detection_limit(self, quantity: str = "c", chip: Chip = None) -> float:
        """Return the D.L. in [mol/s] (quantity='n_dot'), [Pa] ('p'), or [mM] ('c')

        If quantity is 'p' or 'c', then it needs a chip with the right T, p, and carrier
        """
        if self.background_std is None:
            return None
        n_dot_lim = 3 * self.background_std / self.F
        if quantity == "n_dot":
            return n_dot_lim
        if not chip:
            raise TypeError(
                "need a chip to calculate detection limit as pressure or concentration"
            )
        T = chip.T
        p = chip.p
        n_dot_0 = chip.calc_n_dot_0()
        if quantity == "p":
            return p * n_dot_lim / n_dot_0
        if quantity == "c":
            H = MoleculeDict().get(self.mol).calc_H(n_dot_0=n_dot_0, T=T, p=p)
            return n_dot_lim / H
        raise NotImplementedError(
            "Options for detection limit quantity are 'n_dot', 'p', and 'c'"
        )


class Calibration(SensitivityList):
    """Class for saving and loading measured sensitivity factors

    Dev time usage is to directly initiate a calibration with a list of CalPoints
    (cal_list argument to __init__()) and then save it (save())
    Quant time usage is to load a calibration (load()) and then use it to generate
    one or more sensitivity matrices (make_sensitivity_matrix())
    """

    # ---- methods whose primary purpose is interface with the .json ---- #

    def __init__(
        self,
        # calibration basics:
        cal_list: Optional[List[CalPoint]] = None,
        name: Optional[str] = None,
        *,
        setup: Optional[str] = None,
        date: Optional[str] = None,
        description: Optional[str] = None,
        default_setting: Optional[str] = None,
        mol_props: Optional[
            Dict[str, Any]
        ] = None,  # will override stored molecule properties
        fit_specs: Optional[Dict[str, float]] = None,
        mdict: Optional[MoleculeDict] = None,
        **kwargs: Any,  # extra stuff, such as external settings during calibration
    ):
        """Create a Calibration object given its properties

        Should only be called directly for a brand-new calibration, with a cal_list as
        the main data-containing input.

        To load an existing calibration, use `Calibration.load`.

        Args:
            cal_list (list of CalPoint): The list of calibrations points (takes the
                place of ``SensitivityList.sf_list``), each having as a minimum a mol, mass
                and sensitivity factor in [C/mol].
            name (str): The name of the calibration. By default,
                name = f"{date}_{setup}_{detector}_{description}_{TODAY}"
            setup (str): The setup for which the calibration was made / is valid
            date (str): The date of the calibration measurements
            description (str): Extra info on calibration conditions / intentions
            default_setting (str): The setting to which `Calibration.fit` will correspond.
                Defaults to the setting of the first cal_point
            mol_props (dict): Optional. Properties (e.g. Hcp, spectrum) to override
                props in the molecule data files. form is ``mol: {prop:val}``
            fit_specs (dict): Mapping of settings to fit spect i.e:
                ``setting: {"alpha": ..., "beta": ...}`` or merely:
                ``{"alpha": ..., "beta": ...}`` in which it will be assigned to the default
                setting
            mdict (MoleculeDict): In dev-time, or if in doubt, don't use this.
                Calibration will initiate the mdict. But in quant time, Quantifier
                initiates the mdict and Calibration updates it.
            kwargs: Additional key-word arguments are stored as self.extra_stuff
                this includes, for example, dev-time external conditions.
        """
        # ------ calibration results ------- #
        if cal_list is None:
            cal_list = []

        # cal_list becomes self.sensitivity_list
        super().__init__(cal_list)  # type: ignore

        # ----- calibration basics -------
        if name is None:
            name = f"{date}_{setup}_{description}_{TODAY}"
        self.name = name
        self.setup = setup
        self.date = date or TODAY
        self.description = description
        if cal_list:
            self.default_setting = default_setting or self.cal_list[0].setting
        else:
            print("Warning!!! Calibration empty.")
            self.default_setting = None

        # ------- the molecules -------- #
        if mol_props is None:
            mol_props = {}
        self.mol_props = mol_props
        self.mdict = self.initiate_mdict(mdict=mdict)  # THE collection of molecules

        # ----- fit, i.e. F vs f stuff -------------- #
        self._fits: Dict[str, SensitivityFit] = {}
        if fit_specs and set(fit_specs.keys()) == {"alpha", "beta"}:
            self.fit_specs = {self.default_setting: fit_specs}
        else:
            self.fit_specs = cast(Dict[str, Dict[str, float]], fit_specs)

        # -------- extra stuff -------- #
        self.extra_stuff = kwargs

    def as_dict(self) -> Dict[str, Any]:
        """Return a dictionary including everything needed to recreate self"""
        cal_dicts = [cal.as_dict() for cal in self.cal_list]
        # calibration basics:
        self_as_dict = cast(
            Dict[str, Any],
            dict(
                name=self.name,
                setup=self.setup,
                date=self.date,
                description=self.description,
            ),
        )
        # initial calibration results:
        self_as_dict.update(mol_props=self.mol_props, cal_dicts=cal_dicts)
        self_as_dict.update(fit_specs=self.fit_specs)
        self_as_dict.update(self.extra_stuff)
        return self_as_dict

    def save(self, file_name: str = None, cal_dir: PATHLIKE = None) -> None:
        """Save the `as_dict` form of the calibration to a .json file

        Args:
            file_name (str): Name of file to save in, must end in .json
                If not specified, will use self.name + ".json"
            cal_dir: Path to directory to save calibration in, defaults to
                :attr:`Config.calibration_directory`
        Raises:
            ValueError: if name not specified and "None" in self.name
                (this is the case if name wasn't specified
                and either date or setup wasn't specified when initiating)
        """
        if file_name is None:
            if "None" in self.name:
                raise ValueError("Calibration.save demands a file_name!")
            else:
                file_name = self.name + ".json"
        cal_dir = cal_dir or CONFIG.calibration_directory
        path_to_file = (Path(cal_dir) / file_name).with_suffix(".json")
        self_as_dict = self.as_dict()
        with open(path_to_file, "w") as json_file:
            json.dump(self_as_dict, json_file, indent=4)

    @classmethod
    def load(
        cls, file_name: PATHLIKE, cal_dir: PATHLIKE = None, **kwargs: Any
    ) -> "Calibration":
        """Loads a calibration object from a .json file

        Args:
            file_name: Name of the calibration file
            cal_dir: Path to directory to save calibration in, defaults to
                :attr:`Config.calibration_directory`
            kwargs: (Other) key word arguments are fed to `Calibration.__init__`

        Returns:
            Calibration: A calibration object ready to quantify your data!
        """
        cal_dir = cal_dir or CONFIG.calibration_directory
        path_to_file = (Path(cal_dir) / file_name).with_suffix(".json")
        with open(path_to_file, "r") as json_file:
            self_as_dict = json.load(json_file)

        cal_dicts = self_as_dict.pop("cal_dicts")

        cal_list = [CalPoint(**cal_dict) for cal_dict in cal_dicts]
        self_as_dict["cal_list"] = cal_list

        if "fit_specs" not in self_as_dict and (
            "alpha" in self_as_dict and "beta" in self_as_dict
        ):
            # this will be the case for old calibration files
            fit_specs = dict(
                alpha=self_as_dict.pop("alpha"), beta=self_as_dict.pop("beta")
            )
            self_as_dict["fit_specs"] = fit_specs

        self_as_dict.update(kwargs)
        return cls(**self_as_dict)

    # ----------- methods to do with adding, filtering, and fitting -----------

    @property
    def cal_list(self) -> List[SensitivityFactor]:
        """``cal_list`` is a pseudonymn for ``sf_list`` (parent class jargon)"""
        return self.sf_list

    @cal_list.setter
    def cal_list(self, cal_list: List[SensitivityFactor]) -> None:
        """``cal_list`` is a pseudonymn for ``sf_list`` (parent class jargon)"""
        self.sf_list = cal_list

    @property
    def mol_list(self) -> MOLLIST:
        """The mol list of a calibration is all molecules for which it has a `CalPoint`"""
        return list({cal.mol for cal in self})

    @property
    def mass_list(self) -> MASSLIST:
        """The mass list of a calibration is all masses for which it has a `CalPoint`"""
        return list({cal.mass for cal in self})

    # Note, __add__ and filter are in SensitivityList, but need to be modified to
    # include more than the cal_list. SensitivityList.append() and __iadd__ are okay
    # as is.

    def __add__(self, other: "Calibration") -> "Calibration":  # type: ignore
        """Adding another calibration just appends to the `cal_list`"""
        cal_list = self.cal_list + other.cal_list
        self_as_dict = self.as_dict()
        self_as_dict.pop("cal_dicts")
        self_as_dict["cal_list"] = cal_list
        return Calibration(**self_as_dict)

    def filter(self, **kwargs: SENSITIVITYLIST_FILTER_TYPE) -> "Calibration":
        """See `SensitivityList.filter`. Calibration's implementation retains metadata"""
        cal_list = SensitivityList.filter(self, **kwargs).sf_list
        self_as_dict = self.as_dict()
        self_as_dict.pop("cal_dicts")
        self_as_dict["cal_list"] = cal_list
        return Calibration(**self_as_dict)

    def get(self, mol: str, mass: str) -> Optional[Union[CalPoint, SensitivityUnion]]:
        """Return the `CalPoint` (1) or `SensitivityUnion` (>1) with cals of mol at mass"""
        cal_list = self.filter(mol=mol, mass=mass).cal_list
        if len(cal_list) == 1:
            return cast(CalPoint, cal_list[0])
        if len(cal_list) > 1:
            sensitivity_union = cal_list[0].union(cal_list[1])
            for cal_point in cal_list[2:]:
                sensitivity_union = sensitivity_union.union(cal_point)
            return sensitivity_union
        return None

    @property
    def setting_list(self) -> List[str]:
        """Return the ``settings`` for each `CalPoint` in this object"""
        return list({cal.setting for cal in self})

    def make_fit(self, setting: str) -> SensitivityFit:
        """Make a `SensitivityFit` for the CalPoints with a particular mass setting"""
        fit_spec = self.fit_specs.get(setting, {}) if self.fit_specs else {}
        fit = SensitivityFit(self, setting=setting, **fit_spec)  # type: ignore
        self._fits[setting] = fit
        return fit

    def make_fits(self) -> None:
        """Make a `SensitivityFit` for the CalPoints for each mass setting"""
        for setting in self.setting_list:
            self.make_fit(setting)

    def get_fit(self, setting: str) -> SensitivityFit:
        """Return the cached or newly made fit at `setting`"""
        if setting in self._fits:
            return self._fits[setting]
        return self.make_fit(setting)

    @property
    def default_fit(self) -> SensitivityFit:
        """The fit at the default setting of the calibration"""
        return self.get_fit(self.default_setting)

    @property
    def fit(
        calibration,
    ) -> Union["_MyMultiSettingFit", SensitivityFit]:
        """Return an object whose functions choose which fit to use based on setting"""
        setting_list = calibration.setting_list
        if len(setting_list) == 1:
            return calibration.get_fit(setting_list[0])

        return _MyMultiSettingFit(calibration)

    # ---- methods whose primary purpose is interface with quant.Molecule ---- #

    def __getitem__(self, key: MOL) -> Molecule:  # type: ignore
        """Indexing a `Calibration` looks up a molecule in ``calibration.mdict``"""
        try:
            return cast(Molecule, self.mdict[key])
        except KeyError:
            raise KeyError(
                f"self.mdict has no key '{key}'. Try self.molecule('{key}')."
            )

    def molecule(self, mol: MOL) -> Molecule:
        """Return a `Molecule` instance, calibrated if possible

        dev-time only

        Returns:
            Molecule: An `Molecule` instance of mol loaded with the available
                and/or requested calibration and quantification values
        """
        if mol in self.mdict:
            return cast(Molecule, self.mdict[mol])
        else:
            return self.mdict.get(mol)

    def initiate_mdict(self, mdict: Optional[MoleculeDict] = None) -> MoleculeDict:
        """Populate ``self.mdict`` based on ``self.cal_list``

        This function is how `Calibration.__init__` loads its molecules and adds
        saved information into them. It can also be used to add more molecules or
        more information into the calibration and its molecules afterwards.

        It uses the singleton class `MoleculeDict` initialize the molecule collection.
        As calibration is typically the first gateway to new or saved sensitivity
        data, this is in dev-time the one and only place in quant that an instance of
        MoleculeDict is created. In quant time, though, `Quantifier` will initiate a
        medium-aware mdict and `Calibration` will update it here.

        Effects:

         * Sets self.mdict to the MoleculeDict

        Returns:
            dict: mdict, which has molecule names as keys and `Molecule` objects as values

        """
        if not mdict:
            mdict = MoleculeDict()
        # ----- update molecules with info in self.cal_list -----
        for cal in self:
            m = mdict.get(cal.mol)
            m.update(F={cal.mass: cal.F})

        # ----- update molecules with info in self.mol_props -----
        for mol, props in self.mol_props.items():
            m = mdict.get(mol)
            for key, value in props.items():
                setattr(m, key, value)

        # Return it rather than set it so that self.mdict is explicitly set in __init__:
        return mdict

    # ------ methods having to do with the sensitivity matrix ------- #

    def make_sensitivity_matrix(
        self,
        mol_list: Optional[MOLLIST] = None,
        mass_list: Optional[MASSLIST] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> SensitivityMatrix:
        """Make the sensitivity matrix from the original calibration points.

        This function probably only gets called on initiation. It also gets called if
        `mol_list` or `mass_list` are changed (shouldn't happen during quant time).

        Args:
            mol_list (list of str): Names of molecules to correspond to rows in sm.F_mat
            mass_list (list of str): The masses to correspond to columns in sm.F_mat
            metadata (dict): Metadata for the sensitivity matrix. By default, it will be
                populated with keys "name", "sm_name", "time", and "number"

        Effects:

         * Sets self.sm to the new sensitivity matrix
         * Adds the new sensitivity matrix to self.sm_dict

        Returns:
            SensitivityMatrix: The new active sensitivity matrix
        """
        metadata_0 = dict(  # wow this is annoying. Can't wait for dict union with "|"
            time=time.time(),
        )
        if metadata:
            metadata_0.update(metadata)
        metadata = metadata_0
        setting_list = list({mass_to_setting(mass) for mass in mass_list})
        if len(setting_list) == 1:
            fit = self.get_fit(setting_list[0])
        else:
            fit = cast(SensitivityFit, self.fit)  # this is a MultiSettingFit imposter
        sm = self.to_sensitivity_matrix(
            mol_list=mol_list, mass_list=mass_list, metadata=metadata, fit=fit
        )
        return sm

    @classmethod
    def load_as_sm(
        cls,
        *,
        file_name: PATHLIKE,
        mol_list: MOLLIST,
        mass_list: MASSLIST,
        cal_dir: PATHLIKE = None,
        **kwargs: Any,
    ) -> SensitivityMatrix:
        """Load the calibration and immediately make a sensitivity matrix with it

        Args:
            file_name (str or Path): Name of the calibration file.
            mol_list (list of str): Names of molecules to correspond to rows in sm.F_mat
            mass_list (list of str): The masses to correspond to columns in sm.F_mat
            cal_dir: Path to directory to save calibration in, defaults to
                :attr:`Config.calibration_directory`
            kwargs: (Other) key word arguments are fed to `Calibration.__init__`

        """
        cal_dir = cal_dir or CONFIG.calibration_directory
        calibration = cls.load(file_name=file_name, cal_dir=cal_dir, **kwargs)
        return calibration.make_sensitivity_matrix(mol_list, mass_list)

    # ------ methods for visualizing and sanity-checking the calibration ------- #

    def fit_F_vs_f(self, setting: str = None) -> None:
        """Shortcut to the fit_F_vs_f of the calibration's `SensitivityFit` at ``setting``"""
        fit = self.get_fit(setting=setting or self.default_setting)
        return fit.fit()

    def plot_F_vs_f(self, setting: str = None, **kwargs: Any) -> "pyplot.Axes":
        """Shortcut to the plot_F_vs_f of the calibration's `SensitivityFit` at ``setting``"""
        fit = self.get_fit(setting=setting or self.default_setting)
        return fit.plot_F_vs_f(**kwargs)

    def fit_all(self) -> None:
        """Fit F vs f for each setting in the calibration"""
        for setting in self.setting_list:
            self.fit_F_vs_f(setting=setting)

    def plot_all(self, **kwargs: Any) -> Dict[str, "pyplot.Axes"]:
        """Plot F vs f for each setting in the calibration"""
        axes = {}
        for setting in self.setting_list:
            ax = self.plot_F_vs_f(setting=setting, **kwargs)
            ax.set_title(setting)
            axes[setting] = ax
        return axes

    def plot_as_spectrum(
        self,
        mol: Optional[MOL] = None,
        ax: Union[str, "pyplot.Axes"] = "new",
        color: Optional[str] = None,
        **kwargs: Any,
    ) -> "pyplot.Axes":
        """Plot each cal_point in self on an F vs m/z plot.

        This works well with filter. For example,
        ``calibration.filter(mol="CH4").filter(setting="CEM").plot_as_spectrum()``
        plots the measured absolute-sensitivity CEM-detector spectrum of CH4.

        """
        if ax == "new":
            fig, ax = make_axis()
            ax.set_xlabel("m/z in atomic units")
            ax.set_ylabel("sensitivity / [C/mol]")
        ax = cast("pyplot.Axes", ax)

        if not mol:
            for mol_ in self.mol_list:
                self.plot_as_spectrum(mol=mol_, ax=ax, color=color, **kwargs)
        else:
            color_i = color or STANDARD_MOL_COLORS.get(mol, None)
            for n, cal_point in enumerate(self.filter(mol=mol)):
                M = mass_to_M(cal_point.mass)
                F = cal_point.F
                ax.plot([M, M], [F, 0], ":", color=color_i)  # type: ignore
                if n == 0:  # only label one so that the axis looks good.
                    ax.plot(M, F, "o", color=color_i, label=mol, **kwargs)  # type: ignore
                else:
                    ax.plot(M, F, "o", color=color_i, **kwargs)  # type: ignore
            ax.legend(loc="upper right")
        return ax

    def prints_report(self, chip: Optional[Chip] = None, long: bool = True) -> str:
        """Return a string which is a report of the calibration with accuracies, etc.

        The precision, accuracy, and detection limit are all determined from the
        single mass with the highest sensitivity factor. If the CalPoint is a
        `SensitivityUnion`, then the accuracy is the relative standard deviation
        of the contained sensitivity factors, and the precision and detection limit
        come from the first contained `CalPoint` with a stored precision.
        The precision must be stored in the `CalPoint` when it is derived from
        experimental data. It is the relative error when the calculated flux and fitted
        sensitivity factor to explain the measured signal.

        The detection limit is determined by the background noise, which must be
        determined during the calibration experiment and stored in the CalPoint. This
        noise is converted first into a flux by dividing by the signal and then into
        a concentration using the mass transfer coefficient determined by the p, T,
        and total flux of a chip.

        Precision, accuracy, and detection limit all implicitly assume there are no
        interferences that were not present during the calibration experiment.

        Args:
            chip (Chip): The chip with ``p``, ``T`` and total flux used to determine D.L.
            long (bool): Whether to include all the individual CalPoints (True)

        Returns:
            str: A report as a single string containing linebreaks with general
            calibration metadata at the top and then a section for each molecule in the
            calibration with the averaged sensitivities, the accuracy, precision, and
            detection limit if available. This is followed by the representation of
            each of the CalPoints available for the molecule.

        """
        sf_dict = cast(Dict[MOL, CalPoint], self.to_sf_dict())
        report_lines = []
        report_lines += [
            f"calibration name = {self.name}\n",
            f"setup = {self.setup}\n",
            f"date = {self.description}\n",
            f"description = {self.description}\n",
            f"includes the following molecules: {self.mol_list}\n",
            "\n",
        ]
        for mol, sf_dict_i in sf_dict.items():
            spectrum = {mass: sf.F for mass, sf in sf_dict_i.items()}  # type: ignore
            mass_max = max([(F, mass) for mass, F in spectrum.items()])[1]
            sf_max = sf_dict_i[mass_max]  # type: ignore
            if isinstance(sf_max, SensitivityUnion):
                accuracy = sf_max.accuracy * 100
                try:
                    sf_max = [sf for sf in sf_max.sf_list if sf.precision][0]  # type: ignore
                except IndexError:
                    sf_max = sf_max.sf_list[0]
            else:
                accuracy = None
            precision = sf_max.precision and sf_max.precision * 100
            if chip and sf_max.background_std:
                detection_limit = sf_max.calc_detection_limit(chip=chip)
                dl_mgpl = detection_limit * self.mdict.get(mol).M
            else:
                detection_limit = None
                dl_mgpl = None

            report_lines += [
                "-" * 100 + "\n",
                mol + "\n",
                f"\tabsolute sensitivities in [C/mol] = {spectrum}\n",
                f"\tSpecs (smaller is better!) are based on {mass_max}:\n"
                f"\t\t(1 - accuracy) in [%] = {accuracy}\n",
                f"\t\t(1 - precision) in [%] = {precision}\n",
                f"\t\tdetection limit in [mM] = {detection_limit}\n",
                f"\t\tdetection limit in [mg/l] = {dl_mgpl}\n",
                "\n",
            ]
            if long:
                report_lines += ["\tAll measurements:\n"]
                for mass, sf in sf_dict_i.items():  # type: ignore
                    if isinstance(sf, SensitivityUnion):
                        for sf_n in sf.sf_list:
                            report_lines += [f"\t\t{sf_n}\n"]
                    else:
                        report_lines += [f"\t\t{sf}\n"]
                report_lines += ["\n"]

        return "".join(report_lines)

    def print_report(self, *args: Any, **kwargs: Any) -> None:
        """Print the report generated by `prints_report` to the terminal"""
        print(self.prints_report(*args, **kwargs))


class _MyMultiSettingFit:
    def __init__(self, calibration: Calibration) -> None:
        self.calibration = calibration

    def predict_sf(self, mol: MOL, mass: MASS) -> SensitivityFactor:
        setting = mass_to_setting(mass)
        # print(f"Using {self}! Predicting {mol} at {mass}, setting={setting}")
        return self.calibration.get_fit(setting).predict_sf(mol, mass)

    def predict_F(self, mol: MOL, mass: MASS) -> float:
        setting = mass_to_setting(mass)
        return self.calibration.get_fit(setting).predict_F(mol, mass)

    def __getattr__(self, item: str) -> Any:
        return getattr(self.calibration.default_fit, item)
