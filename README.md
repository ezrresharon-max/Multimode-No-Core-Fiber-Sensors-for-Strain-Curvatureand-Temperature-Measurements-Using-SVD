# Calibration and Demodulation Instructions

This repository provides experimental spectral data from three fabricated multimode–no-core–multimode (MNM) fiber-optic sensors, designated MNM-1, MNM-2, and MNM-3. For each sensor, the repository includes spectra measured in response to changes in curvature, strain, and temperature, along with MATLAB scripts for calibration, demodulation, and cross-sensitivity analysis.

## Repository Structure

Repository root/
├── 1. README.txt
├── 2. MATLAB_CODE/
│   ├── crosssens.m
│   ├── demodula.m
│   └── svdsystem.m
└── 3. EXPERIMENTAL_DATA/
    ├── CURVA960/
    ├── CURVA965/
    ├── CURVA977/
    ├── STRAIN960/
    ├── STRAIN965/
    ├── STRAIN977/
    ├── TEMP960/
    ├── TEMP965/
    └── TEMP977/

| Sensor | Dataset suffix | Curvature | Strain | Temperature |
|---|---|---|---|---|
| MNM-1 | 960 | `CURVA960` | `STRAIN960` | `TEMP960` |
| MNM-2 | 965 | `CURVA965` | `STRAIN965` | `TEMP965` |
| MNM-3 | 977 | `CURVA977` | `STRAIN977` | `TEMP977` |

Each dataset folder contains the experimental spectra for one sensor and one measurand. The MNM-3 scripts expect five measurement rounds in subfolders named `1ra medicion`, `2da medicion`, `3ra medicion`, `4ta medicion`, and `5ta medicion`.

## MATLAB Setup

Set MATLAB's **Current Folder** to the repository root and add the scripts to the MATLAB search path:

```matlab
addpath(fullfile(pwd, '2. MATLAB_CODE'))
```

The supplied `svdsystem.m` and `crosssens.m` files still use relative paths from the previous folder structure. Update the assignments shown below before running those scripts.

## SVD-Based Calibration (`svdsystem.m`)

The supplied configuration uses the MNM-3 datasets. In the configuration block of `svdsystem.m`, update the three `dataFolder` assignments:

```matlab
ana(1).dataFolder = fullfile('3. EXPERIMENTAL_DATA', 'CURVA977');
ana(2).dataFolder = fullfile('3. EXPERIMENTAL_DATA', 'STRAIN977');
ana(3).dataFolder = fullfile('3. EXPERIMENTAL_DATA', 'TEMP977');
```

With the repository root still set as the **Current Folder**, run `svdsystem` from the MATLAB Command Window. The script extracts the applied values from the spectral filenames and runs without interactive input.

The script uses rounds T1–T3 to construct a rank-four singular value decomposition (SVD) basis and sixth-order polynomial calibration trajectories over 840–1000 nm. It uses rounds T4 and T5 for validation.

The figures show the experimental spectra, SVD components and amplitudes, normalized coefficient fits, three-dimensional calibration trajectories, and estimated versus applied values. The summary reports the validation root mean square error (RMSE), along with the coefficient of determination (R²) and standard deviation from the calibration regression.

The script saves the generated `.mat` files in MATLAB's **Current Folder**, which is the repository root when following the setup above. When running `demodula`, select the calibration files generated in that folder.

- `calibracioncurvatura_977.mat`
- `calibracionstrain_977.mat`
- `calibraciontemperatura_977.mat`

To calibrate MNM-1 or MNM-2, select the corresponding dataset folders from the table and check that the `fileGlob` and `regex` settings match the actual spectral filenames. Use distinct `matOut` filenames with suffix `960` or `965` to keep the calibration models for each sensor identifiable.

## Spectral Demodulation (`demodula.m`)

After generating the calibration models, run `demodula` from the MATLAB Command Window. Select the three calibration files for the same sensor in this order:

1. Curvature (`calibracioncurvatura_977.mat` for MNM-3).
2. Strain (`calibracionstrain_977.mat` for MNM-3).
3. Temperature (`calibraciontemperatura_977.mat` for MNM-3).

Then select one or more experimental spectra from that sensor. For each spectrum, the script identifies the measurand by selecting the calibration model with the smallest normalized residual and estimates the corresponding value.

If the script can extract the applied value from the filename, the results table also includes the reference value and estimation error. The script plots the last analyzed spectrum relative to the three calibration trajectories.

## Cross-Sensitivity Analysis (`crosssens.m`)

This script compares the MNM-3 spectral responses to curvature, temperature, and strain in a common SVD basis. In `crosssens.m`, update the `ddir` assignment:

```matlab
ddir = fullfile(pwd, '3. EXPERIMENTAL_DATA');
```

With the repository root set as the **Current Folder**, run `crosssens` from the MATLAB Command Window. The script reads the five measurement rounds directly from `CURVA977`, `TEMP977`, and `STRAIN977` under `3. EXPERIMENTAL_DATA`.
