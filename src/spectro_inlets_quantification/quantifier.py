# This file is under dual PROPRIETARY and GPL-3.0 licenses. See DUAL_LICENSE for details.

"""Outer level run-time module

This is the master module for quantification. It uses the sensitivity factors in the
SensitivityMatrix of a Calibration, as well as the capillary function and gas properties
in Chip and Gas, and the external conditions in Medium, to define a Quantifier whose
method calc_c(signal) turns data processed by SignalProcessor into concentrations.

"""


from .medium import Medium
from .molecule import MoleculeDict
from .chip import Chip
from .calibration import Calibration
from .config import Config


CONFIG = Config()


class Quantifier:
    """Class for quantifying signals to fluxes or concentrations"""

    def __init__(
        self,
        calibration=None,
        calibration_file=None,
        cal_dir=None,
        # quantification basics:
        mol_list=None,  # list of names of molecules to be calibrated
        mass_list=None,  # list of masses for which signals will be provided
        mol_and_mass_list_pairs=None,
        sensitivity_matrices=None,  # dict of sensitivity matrices
        # quantification measurement external conditions:
        T=None,  # temperature during quantification
        p=None,  # pressure during quantification
        carrier=None,  # carrier gas composition during quantification measurement
        chip="standard.json",  # object of quant.Chip or chip file name
        pp_mode="He_solver",
        relaxed=None,  # molecules with flexible flux when calculating partial p,
        verbose=True,
    ):
        """Create a Quantifier object given its properties

        Args:
            calibration (Calibration): the calibration, if already loaded.
            calibration_file (string): the file from which to load the calibration, if
                a pre-loaded calibration is not given
            cal_dir: The directory from which to load calibrations, defaults to
                :attr:`Config.calibration_directory`
            mol_list (list of str):
                names of molecules to be quantified by the calibration.
            mass_list (list of str):
                masses that the calibration will use to quantify.
            mol_and_mass_list_pairs (list): List of (mol_list, mass_list) tuple pairs,
                each constituting a sensitivity matrix
            sensitivity_matrices (dict): A dictionary of sensitivity matrices. If this
                argument is given, these sensitivity matrices will be used instead of
                building new ones from `calibration` and mol and mass lists. Keys in the dict
                are 'sm_list' and 'master_sm'
            carrier: (dict or str): composition of the carrier gas during
                quantification measurements. Possible values the same as for carrier0.
            T (float): The temperature during calibration measurements in [K]
                Defaults to STANDARD_TEMPERATURE in .medium
            p (float): The pressure the calibration measurements in [Pa]
                Defaults to STANDARD_PRESSURE in .medium
            chip (Chip or str or dict): the chip used during quantification
                measurements. Can be an instance of the Chip class or a string
                specifying the file describing the chip or kwargs for Chip.__init__
            pp_mode (str): The way in which chip partial pressures are calculated from
                MS-calibrated flux and evt. the capillary equation. See 'calc_pp'
            relaxed (set of str): names of molecules that can have a different flux
                when calculating the capillary gas composition than implied by
                quantification of the MS signals. By default this will include the
                component(s) of the carrier gas. See Chip.partial_pressures_by_mix_in
                for details.
            verbose (bool): Whether to make output verbose
        """
        # ------- quantification external conditions ------ #
        self.medium = Medium()  # initiates the one-and-only medium
        # FIXME these looks like bugs, they should save T and p on medium, probably shows that
        # they have never been used
        if T:
            self.T = T
        if p:
            self.p = p
        if isinstance(chip, str):
            # interpret as a file name, with the .json extension!
            chip = Chip.load(chip, verbose=False)
        elif isinstance(chip, dict):
            # then interpret it as a Chip.as_dict() dictionary (same as json contents)
            chip = Chip(**chip)
        self.chip = chip
        self.chip.carrier = carrier  # this will make sure it's a Gas object.
        chip.medium = self.medium
        self.mdict = MoleculeDict(medium=self.medium)

        # ------- the calibration and SensitivityMatrix -------------- #
        cal_dir = cal_dir or CONFIG.calibration_directory
        self.calibration_file = calibration_file
        self.calibration = calibration or Calibration.load(
            calibration_file,
            cal_dir=cal_dir,
            mdict=self.mdict,
        )

        if sensitivity_matrices:
            self._sm_list = sensitivity_matrices["sm_list"]
            self._master_sm = sensitivity_matrices["master_sm"]
        else:
            self._sm_list = []
            if mol_list and mass_list:
                sm = self.make_sm(mol_list=mol_list, mass_list=mass_list)
                self._sm_list += [sm]
            if mol_and_mass_list_pairs:
                for mol_list_n, mass_list_n in mol_and_mass_list_pairs:
                    self._sm_list += [self.make_sm(mol_list_n, mass_list_n)]
            self._master_sm = None

        # signal and concentration are set and calculated, respectively, later:
        self.S = {}  # signals / [A]

        self.pp_mode = pp_mode

        if relaxed is None:
            relaxed = self.chip.carrier.mol_list
        self.relaxed = set(relaxed)

        self.verbose = verbose

        # comp and n_dot and mass transfer coefficients are calculated later:
        self.n_dot = {}
        self.n_dot_0 = None
        self.n_dot_analyte = {}
        self.H = {}
        self.p_app = None

    def __eq__(self, other):
        """Return whether a quantifier is equal"""
        if type(self) != type(other):
            return False
        if len(self.sm_list) != len(other.sm_list):
            return False
        if self.chip != other.chip:
            return False
        return self.sm_list == other.sm_list

    @classmethod
    def from_dict(cls, obj_as_dict, chip, carrier, calibration_file):
        """Create object from dictionary containing the object state

        Args:
            object_as_dict (dict): The dict representation of the object, as returned by
                `as_dict`

        Remaining arguments as described in `__init__`.

        TODO: This mirrors the half implementation in to_dict and so requires extra
        mandatory arguments

        TODO: The dict is, at this point, not JSON able as it should be, but it is pickable.

        """
        if "sensitivity_matrices" in obj_as_dict:
            from .sensitivity import SensitivityMatrix

            sensitivity_matrices = obj_as_dict["sensitivity_matrices"]
            sensitivity_matrices["sm_list"] = [
                SensitivityMatrix(**sensitivity_matrix_def)
                for sensitivity_matrix_def in sensitivity_matrices["sm_list"]
            ]
            sensitivity_matrices["master_sm"] = SensitivityMatrix(
                **sensitivity_matrices["master_sm"]
            )
        return cls(
            **obj_as_dict,
            chip=chip,
            carrier=carrier,
            calibration_file=calibration_file,
        )

    def to_dict(self):
        """Return the object state as a dictionary

        TODO: This is only a partial implementation, saving the most pertinent information
        items, but leaving a few out which will need to be re-supplied when re-creating the
        object with :meth:`from_dict`

        """
        return {
            "sensitivity_matrices": {
                "sm_list": [
                    sensitivity_matrix.as_dict() for sensitivity_matrix in self.sm_list
                ],
                "master_sm": self.master_sm.as_dict(),
            },
            "pp_mode": self.pp_mode,
            "verbose": self.verbose,
        }

    def make_sm(self, mol_list, mass_list):
        """Shortcut to the Calibration's make_sensitivity_matrix"""
        return self.calibration.make_sensitivity_matrix(
            mol_list=mol_list, mass_list=mass_list
        )

    @property
    def sm(self):
        """The primary SensitivityMatrix"""
        if self._sm_list:
            return self._sm_list[0]

    @property
    def sm_list(self):
        """list of SensitivityMatrix"""
        return self._sm_list

    @property
    def master_mass_list(self):
        """Union of the mass lists of all the quantifier's sensitivity matrices"""
        mass_list = set()
        for sm in self.sm_list:
            for mass in sm.mass_list:
                mass_list.add(mass)
        return list(mass_list)

    @property
    def master_mol_list(self):
        """Union of the mol lists of all the quantifier's sensitivity matrices"""
        mol_list = set()
        for sm in self.sm_list:
            for mol in sm.mol_list:
                mol_list.add(mol)
        return list(mol_list)

    @property
    def master_sm(self):
        """Union of all the quantifier's sensitivity matrices"""
        if not self._master_sm:
            self._master_sm = self.calibration.make_sensitivity_matrix(
                mass_list=self.master_mass_list, mol_list=self.master_mol_list
            )
        return self._master_sm

    @property
    def mol_list(self):
        return self.sm.mol_list

    @property
    def mass_list(self):
        return self.sm.mass_list

    @property
    def carrier(self):
        return self.chip.carrier

    @property
    def T(self):
        return self.medium.T

    @T.setter
    def T(self, T):
        self.medium.T = T

    @property
    def p(self):
        return self.medium.p

    @p.setter
    def p(self, p):
        self.medium.p = p

    @property
    def p_vac(self):
        return self.medium.p_vac

    @p_vac.setter
    def p_vac(self, p_vac):
        self.medium.p_vac = p_vac

    # ---- methods for calculating flux ------ #

    def print_sensitivity_matrices(self):
        print(self.prints_sensitivity_matrices())

    def prints_sensitivity_matrices(self):
        """Format all information about the sensitivity matrices and return as a string"""
        out = "Quantifier got these sensitivity matrices:"
        for i, sm in enumerate(self.sm_list):
            out += f"\nsm # {i}:"
            out += f"\n{sm}"
            out += "\n" + sm.prints_F_mat()
        return out

    # ------- methods for calculating n_dot, p, and c ------------ #

    def calc_n_dot(self, signals, sequence=None):
        """Return {mol: n_dot_i} of flux in [mol/s] given {mass: S_M} of signal in [A].

        This is where the multiple sensitivity matrices are used sequentially.

        Each application of the sensitivity matrix will ignore the flux of molecules
        for which the flux was calculated by a previous sensitivity matrix. The
        signal explained by previously-calculated fluxes will be subtracted before
        application of the sensitivity matrix.
        TODO: It might be nice to make it an option to count the flux from the
            "leftover signal" towards the returned flux. That might complicate the
            implementation, so I think we should get a better feel for how it works
            as is first.

        Args:
            signals (dict or SignalDict): {mass: S_M}, where S_M is the signal in [A]
            sequence (list of int): The indeces specifying the order in which to apply
                the sensitivity matrices in self.sm_list. Defaults to the order that
                they appear in self.sm_list.
        """
        n_dot = {}

        sequence = sequence or range(len(self.sm_list))

        for n in sequence:
            sm = self.sm_list[n]
            if self.verbose:
                print(
                    f"\n### Quantifier.calc_n_dot is applying the #{n} " f"sm: {sm} ###"
                )
            all_explained_signals = self.master_sm.calc_signal(n_dot=n_dot)
            explained_signals = {
                mass: all_explained_signals[mass] for mass in sm.mass_list
            }
            unexplained_signals = {
                mass: signals[mass] - explained_signals[mass] for mass in sm.mass_list
            }
            if self.verbose:
                print(f"Signal in [A] already explained = {explained_signals}")
                print(f"Signal in [A] left to explain = {unexplained_signals}")
            n_dot_n = sm.calc_n_dot(signals=unexplained_signals)

            for mol, n_dot_i in n_dot_n.items():
                if mol in n_dot:
                    if self.verbose:
                        print(
                            f"already got n_dot['{mol}'] = {n_dot[mol]*1e9} [nmol/s]"
                            f"but it looks from the unexplained signal "
                            f"like there's an extra {n_dot_i*1e9} [nmol/s]"
                        )
                else:
                    n_dot[mol] = n_dot_i
                    if self.verbose:
                        print(f"got n_dot['{mol}'] = {n_dot_i*1e9} [nmol/s]")
        return n_dot

    def calc_pp(self, signals=None, p=None, T=None, sequence=None, mode=None):
        """calculate the concentration [mol/m^3] for each quantified molecule

        Args:
            signals (dict or SignalDict): {mass: S_M}, where S_M is the signal in [A]
            p (float): pressure [Pa] if to be updated
            T (float): temperature [K] if to be updated
            sequence (list of int): The indeces of the matrices to use. See calc_n_dot.
            mode (str): how to convert flux to partial pressure. See Chip.calc_pp.

        Returns:
            dict: {i: c^i} where c^i [mol/m^3] is the concentration for each
                molecule i in self.mol_list
        """
        if p is not None:
            self.p = p
        if T is not None:
            self.T = T

        n_dot = self.calc_n_dot(signals=signals, sequence=sequence)  # fluxes in [mol/s]

        mode = mode or self.pp_mode
        pp = self.chip.calc_pp(n_dot=n_dot, mode=mode)  # partial pressures in [Pa]

        return pp

    def calc_c(self, signals=None, p=None, T=None, sequence=None, mode=None):
        """calculate the concentration [mol/m^3] for each quantified molecule

        Args:
            signals (dict or SignalDict): {mass: S_M}, where S_M is the signal in [A]
            p (float): pressure [Pa] if to be updated
            T (float): temperature [K] if to be updated
            sequence (list of int): The indeces of the matrices to use. See calc_n_dot.
            mode (str): how to convert flux to partial pressure. See Chip.calc_pp.

        Returns:
            dict: {i: c^i} where c^i [mol/m^3] is the concentration for each
                molecule i in self.mol_list
        """
        pp = self.calc_pp(
            signals=signals, p=p, T=T, sequence=sequence, mode=mode
        )  # partial pressures in [Pa]

        c = {}
        for mol, pp_i in pp.items():
            c[mol] = pp_i / self.mdict[mol].calc_KH()  # concentration in [mol/m^3]

        return c
