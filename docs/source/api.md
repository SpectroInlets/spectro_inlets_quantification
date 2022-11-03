(API)=
# API

```{toctree}
---
maxdepth: 2
caption: Contents
glob:
---
api/*
```

(sec:modules)=
## Overview of the modules of [spectro_inlets_quantification](API)

The [spectro_inlets_quantification](API) package tries to make the best of modular, object-oriented programming to capture all the physics involved in quantitative mass spectrometry with a Spectro Inlets inlet system. Modules and classes do their best to represent physically intuitive concepts in the system and in your workflow. This section gives a very brief description of each module as it exists in quant.physics v1.0, with examples of some of the most important interfaces in the accompanying Ipython notebook.

Most of the equations in the previous sections of this document appear somewhere in the package. It is much bigger than that, though, because of the need to keep various pieces of information in the right places accessible through intuitive interfaces.

The code is pretty well self-documented, so don't hesitate to import anything and input it to python's `help` function to see the docstring. There are code examples to go along with these explanations for each module in the accompanying Ipython notebook. The modules have an explicit hierarchy, and this section goes through them from lowest to highest in that hierarchy.

(sec:config)=
### The {mod}`.config` module

This module contains the singleton {class}`.Config` class, which contains all package wide configurable setting. Importantly, this allows the configuration of the base **data directory**, from which all other data directories (e.g. {attr}`.Config.molecule_directory`, which contains the data files that define e.g. the Henry's-Law parameters for each molecule) are derived.

(sec:constants)=
### The {mod}`.constants` module

The **{mod}`.constants`** module defines constants used throughout [spectro_inlets_quantification](API), such as:

* The physical constants used in [spectro_inlets_quantification](API) such as the `GAS_CONSTANT` representing $R=8.3143$ [J/mol/K].
* The chip design parameters including `STANDARD_CAPILLARY_WIDTH = 6e-6` [m]
* Standard conditions including `STANDARD_TEMPERATURE = 298.15` [K]
* Plotting preferences such as `STANDARD_COLORS`

(sec:tools)=
### The {mod}`.tools` module

The **tools** module contains useful pythony stuff used in multiple places like a {class}`.Singleton` metaclass.

It also includes some little function for parsing mass strings. For example:

```python
>>> from spectro_inlets_quantification.tools import mass_to_M, mass_to_setting
>>> mass_to_M("M32-CEM")
32.0
>>> mass_to_setting("M32-CEM")
'CEM'
```

It also includes a string ({data}`.tools.TODAY`) made when [spectro_inlets_quantification](API) is imported that represent the day's date in Soren's format. For example, the 13th of November 2020 is "20K13" where 'K' (the 11th letter) represents November (the 11th month). This string is used by default when new calibrations are defined.

(sec:medium)=
### The {mod}`.medium` module

The **medium** module's only purpose is to give a home to the one-and-only pressure and temperature of the system, so that different classes agree on these two external conditions. That home, {class}`.Medium`, is a singleton. Accessing the attributes `T` or `p` of ***any*** other class in any of the modules of quant including {class}`.Quantifier`, {class}`.Molecule`, {class}`.Chip`, or {class}`.Gas` will return {attr}`.Medium.T` or {attr}`.Medium.p`, respectively.

(sec:molecule)=
### The {mod}`.molecule` module

The **molecule** module defines the {class}`.Molecule` class. Instances of this class, instantiated with the {meth}`.Molecule.load` (alternative constructor) class method which takes the molecule name (e.g. "H2"), load a number of constants including Henry's-Law constant, thermochemistry, molecular diameter, ionization cross-section data, and reference spectrum from a .json data file in the {attr}`.Config.molecule_directory`.

The {class}`.Molecule` thus instantiated wraps this data with some useful methods, such as one to calculate the volatility constant as a function of temperature (used to generate the data for {numref}`Figure %s <fig:KH>`) and plot the reference spectrum (such as {numref}`Figure %s <fig:NIST>`.

The module also defines a {class}`.MoleculeDict`, a singleton whose only instance is called `mdict` everywhere it is used in [spectro_inlets_quantification](API). The `mdict` collects {class}`.Molecule` instances when they are instantiated with the `mdict.get()` method which takes the molecule name as its argument, and then makes them available by method or by indexing as well. {meth}`.MoleculeDict.get` is almost always preferred over {meth}`.Molecule.load` because it avoids an unnecessary second read of the molecule file.

(sec:mixture)=
### The {mod}`.mixture` module

The **mixture** module has a base class {class}`.Mixture` which serves as a framework for dealing with a mixture of molecules. The most important attribute of its instances (`mixture`) is {attr}`.Mixture.comp`, which is a dictionary of {$i$: $x^i$} where $i$ is the name of a molecule and $x^i$ is its mol fraction in the mixture. If `abc` is an attribute of {class}`.Molecule`, then `mixture.abc` returns the mol-weighted average of `abc` for each of the molecules in the mixture. {meth}`.Mixture.make` is a constructor that takes a `comp` dictionary or the name of a molecule or the name of a standard mixture (such as `"air"`), and populates the mixture accordingly.

The {class}`.Gas` class inherits from {class}`.Mixture` and adds a couple of things: First, it has an updated viscosity correction that overrides the weighted average inherited from {class}`.Mixture` and instead uses the algorithm in Davidson1993. This was used to generate the data for {numref}`Figure %s <fig:eta>`. Second, {meth}`.Gas.saturated_with` takes a molecule as input and returns a gas that contains that molecule at the mol fraction dictated its vapor pressure and the system pressure.

(sec:chip)=
### The {mod}`.chip` module

The **chip** module defines the {class}`.Chip` class, which is basically a wrapper around the capillary equation, Equation {eq}`eq:capillary`. An instance of {class}`.Chip` can be defined with the {meth}`.Chip.load` method which takes as its argument the name of a .json file in the {attr}`.Config.chip_directory`, which contains the capillary dimensions if they differ from the defaults of a known chip. It can also be initialized with the default dimensions directly from `chip=Chip()`. Either can take `p`, `T`, `carrier`, and `solvent` as arguments, where the first two set the system {attr}`.Medium.p` and {attr}`.Medium.T`, respectively.

{attr}`.Chip.gas` represents the gas in the chip. By default, it is {attr}`.Chip.carrier` saturated with {attr}`.Chip.solvent`. The capillary equation is called with {meth}`.Chip.calc_n_dot_0`. It can take a `gas`, `p`, and `T` as inputs but by default uses those of the chips. This method was used to generate the data in {numref}`Figure %s <fig:cap>`. It also has a {meth}`.Chip.calc_n_dot` function which returns a dictionary giving the flux in [mol/s] for each of the molecules in its gas.

The {class}`.Chip` also has a number of methods for calculating the partial pressure of the gases in it given a set of quantified fluxes, including that which solves Equation {eq}`eq:solver1`. This resets the chip's gas.

(sec:peak)=
### The {mod}`.peak` module

The **peak** module defines a {class}`.Peak` base class and classes that inherit from it (so far only {class}`.GaussPeak`). Each of these classes is initiated with x-and-y data, and have methods for extracting a `height`, `width`, and `center` from it. For the base class, these are extracted by simple means (for example, height is the maximum y value), while for a `gauss_peak` (instance of {class}`.GaussPeak`), they come from the Gaussian fit. A peak that fails to fit raises a {class}`.PeakFitError`. A peaks also has a `signal` attribute, which is the height unless it is flagged with `peak.error=True`, in which case `peak.signal` returns 0.

Each {class}`.Peak` class has a `plot` method for visualization. {numref}`Figure %s <fig:peak>` was made with `gauss_peak.plot()`.

(sec:signal)=
### The {mod}`.signal` module

The **signal** module contains three classes for organizing, visualizing, and analyzing raw data.

The {class}`.SignalDict` class stores signals, as defined in {prf:ref}`def-signal`. When a signal is added, either its time is also added or it gets timestamped, so that {class}`.SignalDict` knows its history. Besides that, it acts like a dictionary which returns the latest signal in \[A\] when indexed with a mass-setting string.

The {class}`.SignalProcessor` class is the main data-processing class of [spectro_inlets_quantification](API). It can be loaded from a processor file containing data on non-linearity and background, or initiated directly. Either way, it can be given a `peak_type` which specifies which {class}`.Peak` class the instance, `signal_processor`, uses. Its {meth}`.SignalProcessor.calc_signal` method class takes in raw x-and-y data, corrects it for non-linearity and/or background, makes a {class}`.Peak`, and then calculates the signal, adds it to its {class}`.SignalDict`, and returns it.

(sec:sensitivity)=
### The {mod}`.sensitivity` module

The **sensitivity** module contains several classes for managing sensitivity factors.

First, the {class}`.SensitivityFactor` class is a wrapper around a single sensitivity factor. A `sensitivity_factor` has the actual number in \[C/mol\] as its attribute {attr}`.SensitivityFactor.F`, with its attributes {attr}`.SensitivityFactor.mol` and {attr}`.SensitivityFactor.mass` specifying the molecule and mass-setting for which it applies. The {class}`.SensitivityUnion` class can unite sensitivity factors with matching mol and mass. Its {attr}`.SensitivityUnion.F` is the average of those of its members, and it has a property {attr}`SensitivityUnion.accuracy` telling the relative deviation thereof.

The {class}`.SensitivityList` is a wrapper around a list of sensitivity factors with a {meth}`.SensitivityList.filter` method that returns a smaller {class}`.SensitivityList` based on any attribute of the `sensitivity_factors`. It also has a {meth}`.SensitivityList.to_sensitivity_matrix` method which takes `mol_list` and `mass_list` as arguments and passes them on, with the needed `sensitivity_factors` to {class}`.SensitivityMatrix`.

The {class}`.SensitivityMatrix` class is the home of the central calculation in quantification, Equation {eq}`eq:Q1` for counting flux, which is {meth}`.SensitivityMatrix.calc_n_dot`. In this method, the inverse to {attr}`.SensitivityMatrix.F_mat` ($\mat{F}$) is taken and matrix multiplied onto the `signals` or `signal_dict` given as an argument, and the result is rearranged as a dictionary `n_dot`. {attr}`SensitivityMatrix.F_mat` is a matrix spanning {attr}`SensitivityMatrix.mol_list` and {attr}`.SensitivityMatrix.mass_list` where any entry that was not available from the `sensitivity_list` that initiated it is predicted by the method described in Section [](sec:f).

That prediction is done by the class {class}`.SensitivityFit`. A `sensitivity_fit` has a {meth}`.SensitivityFit.fit` method which determines the parameters {attr}`.SensitivityFit.alpha` and {attr}`.SensitivityFit.beta` ($\alpha$ and $\beta$ in Equation {eq}`eq:prediction`), and a {meth}`.SensitivityFit.plot_F_vs_f` method for visualization ({numref}`Figure %s <fig:Fvf>`).

(sec:cal)=
### The {mod}`.calibration` module

The **calibration** module defines two classes, each of which inherit from classes in the `sensitivity` module.

{class}`.CalPoint` inherits from {class}`.SensitivityFactor` and adds to its parent metadata attributes and specs such as {attr}`.CalPoint.precision` and {attr}`.CalPoint.background_std` where the latter can be used to calculate its detection limit. A {class}`.CalPoint` should be initialized directly with all of these attributes, which are descriptions of and results from a calibration experiment.

{class}`.Calibration` inherits from {class}`.SensitivityList` and adds to its parent methods for visualization, saving, loading, and fitting. {class}`.Calibration` is the class for collecting and saving calibration results for later use. During quantification, a {class}`.Calibration` is loaded and used to generate `sensitivity_matrices` with an enhanced version of the inherited {meth}`.Calibration.make_sensitivity_matrix` function. {attr}`.Calibration.fit`, whose parameters are saved and loaded with the {class}`.CalPoints`, is passed on to the sensitivity matrices that `calibration` makes. The fit is mainly to be able to sanity-check the `calibration`, usually after `filter`'ing for a specific setting. The {meth}`.Calibration.plot_as_spectrum` method makes figures like {numref}`Figure %s <fig:cal>`. It also has a {meth}`.Calibration.print_report` method, which generates a text output report.

(sec:quantifier)=
### The {mod}`.quantifier` module

Finally, the **quantifier** module defines the {class}`.Quantifier` class, which is initiated with a `calibration_file`, a `mol_list` and `mass_list`, and a `chip`, as well as a few other options. The `quantifier` loads the `calibration` from the file and immediately uses it to build a `sensitivity_matrix` based on the `mol_list` and `mass_list`.

Most importantly, `quantifier` binds the methods for quantification according to the procedure in Section [](sec:solve):

* {meth}`.Quantifier.calc_n_dot` for flux quantification,
* {meth}`.Quantifier.calc_pp` for partial pressure quantification, and
* {meth}`.Quantifier.calc_c` for concentration quantification.

That concludes our tour of [spectro_inlets_quantification](API), most of which is in the Ipython notebook. We hope you enjoyed it and find the tools and physics here useful!
