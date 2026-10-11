# UnidirectionalAAmodel

This repository contains MATLAB code and data related to the study of unidirectional AAmode, focusing on different localization phenomena such as Scale-Free Localization (SFL) and Anderson Localization (AL). The code implements various numerical calculations and visualization techniques under different boundary conditions, including Generalized Boundary Conditions (GBC).

## Repository Structure

### Code Files
- `EcData1.m` - Data processing and computation script.
- `Eigenstates_AL_bd.m` - Calculates eigenstates for Anderson Localization under specific boundary conditions.
- `Eigenstates_SFL_bd.m` - Computes eigenstates for Scale-Free Localization under boundary conditions.
- `Energylevel_gbc.m` - Computes energy levels under Generalized Boundary Conditions.
- `PlotMosaicME.m` - Generates mosaic plots for visualization.
- `Plot_GBC_AL_closevalue.m` - Plots specific eigenvalues close to a given value under GBC and AL.
- `Sfl_AL_GBC_Spectrum.m` - Computes and visualizes the spectrum for SFL and AL under GBC.
- `Typical_Eigenstates_ALvpa1.m` - Extracts and plots typical eigenstates for Anderson Localization.
- `Typical_Eigenstates_SE.m` - Computes and analyzes eigenstates under specific settings.
- `Typical_Eigenstates_Sfl_LE_scaling.m` - Investigates eigenstate scaling for SFL.
- `uniAA_AL_GBC_spectrum_bc.m` - Computes spectrum for AL under GBC with specific boundary conditions.
- `uniAA_AL_GBC_spectrum_gap.m` - Analyzes spectral gaps in AL under GBC.

### Data and Figures
The repository contains multiple folders that store data and figures for visualization:
- `AllDataPlot/` - Contains all generated plots.
- `Fig2c1/`, `Fig2c2/`, `Fig2d/` - Store figure data for specific cases.
- `MosaicData/` - Stores mosaic visualization data.
- `fig2a/`, `fig2b/` - Additional figure datasets.

## Usage
1. Ensure you have MATLAB installed.
2. Run the `.m` scripts to compute eigenstates, spectra, and other relevant properties.
3. Use the plotting scripts to visualize the results.
4. Data files in the respective folders provide input/output for further analysis.

## License
This project is for academic and research purposes. If you use the code or data, please cite appropriately.

## Contact
For questions or collaborations, please contact the repository maintainer.


## Exact paths and compatibility entries

The MATLAB source names above live in `Code/code/`. Figure datasets live
under `Data/Data/` (including `AllDataPlot/` and the individual figure
folders). High-precision prototypes and CSV inputs live in `Others/Others/`.
These historical doubled directory names are retained.

Two identical source pairs now have one maintained implementation each:

| Compatibility filename | Canonical implementation |
| --- | --- |
| `Code/code/Plot_GBC_AL_closevalue.m` | `Code/code/uniAA_AL_GBC_spectrum_bc.m` |
| `Others/Others/Test1.cpp` | `Others/Others/Typical_Eigenstates.cpp` |

The MATLAB alias invokes the canonical script. Relative input paths supplied
by the working directory retain priority. When those inputs are absent, the
default DataAnalysisC mode falls back to the recorded Data1.csv, Data2.csv and
Data3.csv by repository-relative path, so copying them into the code folder
is unnecessary. Alternate analysis modes still require their external data.

```matlab
run(fullfile(pwd,'Code','code','Plot_GBC_AL_closevalue.m'));
```

To use regenerated CSV files in another working directory, add Code/code to
the MATLAB path and invoke the script by name from the data directory. The
compatibility alias preserves that working directory:

```matlab
addpath('/path/to/UnidirectionalAAmodel/Code/code');
cd('/path/to/your/csv-inputs');
Plot_GBC_AL_closevalue;
```

The C++ alias includes the canonical translation unit; compile only one of the
two entry files per executable. It requires MPLAPACK's MPFR backend and the
MPFR C++ interface, neither bundled here. Some MATLAB files using `vpa`
require Symbolic Math Toolbox. The Python prototype imports NumPy, Matplotlib
and mpmath.

The C++ source is a historical prototype: its separate Ev buffer is not
filled by the eigensolver call, and its symmetric solver and boundary-matrix
convention need a dedicated numerical review before its printed eigenvectors
are used. Source deduplication does not validate that prototype. No figures,
datasets, binaries, or historical entry files were deleted.

