# **MOOSE**
**Mock Observation Of Synchrotron Emission**

MOOSE is an interactive tool designed for processing simulated synchrotron emission data. It provides functionalities for calculating intensity, polarization, and Faraday rotation maps, while supporting various physical and instrumental parameters.

---

## **Table of Contents**
1. [Main Features](#main-features)
2. [Context](#context)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Caveats](#caveats)
6. [Contributors](#contributors)

---

## **Main Features**
- Process synchrotron simulations based on mock data.
- Compute different parameters maps like RM or Maximum polarized intensity map.
- Compute Stokes parameters I, Q, U.
- Perform Rotation Measure Synthesis (RM Synthesis) for Faraday rotation studies.
- Handle instrumental parameters such as frequency resolution and simulation box sizes.

---
## **Installation**
### Prerequisites
- Ensure **Julia** (v1.10+) is installed. Download it [here](https://julialang.org/downloads/).
- Compatible with Linux, macOS, and Windows.

### Steps
1. Clone the repository:
    ```bash
    git clone https://github.com/username/MOOSE.git
    cd moose
    ```

2. Install required Julia dependencies:
    ```julia
    using Pkg
    Pkg.add([
        "FITSIO", "DataFrames", "CSV", "Interpolations", "Dierckx",
        "LinearAlgebra", "SpecialFunctions", "QuadGK", "StatsBase",
        "Statistics", "KernelDensity", "StringEncodings", "FFTW",
        "ImageFiltering", "Random", "Distributions"
    ])
    ```
---
## **Context: Faraday Tomography in the ISM**

Faraday tomography is a powerful technique for studying the structure and composition of the interstellar medium (ISM) through its magnetic field and ionized gas. This method leverages the **Faraday rotation effect**, which occurs when a linearly polarized electromagnetic wave passes through a magnetized ionized medium, causing the polarization plane to rotate. The amount of rotation is proportional to the **Faraday depth**, which depends on the strength and direction of the magnetic field, the electron density, and the path length.

### **What Does Faraday Tomography Calculate?**
Faraday tomography reconstructs the polarized emission as a function of Faraday depth (in units of rad m\(^{-2}\)). By performing a **Rotation Measure (RM) synthesis** (Brentjens & de Bruyn, 2005), it disentangles polarized emission from sources at different Faraday depths along the line of sight. The resulting **Faraday dispersion function** maps polarized emission as a function of Faraday depth, providing insights into:
- The magnetic field structure along the line of sight.
- The distribution and density of ionized gas.
- The contribution of multiple layers of emission within the ISM.

For a detailed introduction to Faraday tomography and RM synthesis, see:
- Brentjens & de Bruyn (2005) ([DOI: 10.1051/0004-6361:20052990](https://doi.org/10.1051/0004-6361:20052990))
- Ferrière et al. (2021): ([DOI: 10.1093/mnras/stab1641](https://doi.org/10.1093/mnras/stab1641))

### **Why Use Faraday Tomography?**
Faraday tomography is critical for understanding the magneto-ionic properties of the ISM, as magnetic fields play a central role in:
- Cosmic-ray transport.
- Star formation processes including collapse of molecular clouds or supernovae remnants evolution.
- The dynamics of the ISM, including turbulence and large-scale flows.

### **Using MOOSE for Faraday Tomography**
MOOSE is designed to process data from numerical MHD simulations and perform Faraday tomography on mock synchrotron emission data. To ensure the tool recognizes and processes the data correctly, it requires the following input files **with specific file names**:

1. **Magnetic field cubes**:
   - \( Bx.fits \): Magnetic field component along the x-axis preferably in microG.
   - \( By.fits \): Magnetic field component along the y-axis preferably in microG.
   - \( Bz.fits \): Magnetic field component along the z-axis preferably in microG.
   
2. **Density cube**:
   - \( density.fits \): density preferably in cm\(^{-3}\).

3. **Temperature cube**:
   - \( temperature.fits \): Temperature preferably in Kelvin.

### **Important Note**
For the code to correctly identify and process these files, they **must** be named exactly as follows:
- `Bx.fits`, `By.fits`, `Bz.fits` for the magnetic field components.
- `density.fits` for the density cube.
- `temperature.fits` for the temperature cube.

If your data cubes are named differently, you need to rename them before running the code. Failure to do so will result in errors, as the tool depends on these specific file names to match the required inputs.

For example, if your original files are named `mag_field_x.fits` and `n.fits`, you should rename them:
```bash
mv mag_field_x.fits Bx.fits
mv n.fits density.fits
```
---

## **Usage**
### **MOOSE**
The primary tool in MOOSE is the **`MOOSE()`** function, which guides users interactively through the processing of simulated data. To launch MOOSE, use:
```julia
include("MOOSE.jl")  # Load the main file
MOOSE()              # Start the interactive tool
```
Below is an example session that demonstrates how to use `MOOSE()` step-by-step.

### **Interactive Session Example**
```julia
MOOSE()
```
1. Base directory for simulations:
```julia
"Enter the base directory for simulations (default: /path/to/default/directory):"
/path/to/simulations
```
2. Select Simulations
```julia
"Available simulations:"
[1] /path/to/simulations/Simulation1
[2] /path/to/simulations/Simulation2
[3] /path/to/simulations/Simulation3
"Do you want to process all simulations or choose specific ones? (Enter 'all' or 'choose') (default: all):"
all
```
3. Define Units and Conversion Factors:
```julia
"Is the unit of magnetic field B in μG (microGauss)? (Y/N) (default: N):"
n
"Enter the conversion factor for magnetic field B to μG (microGauss): (default: 1000.0):"

"Is the unit of number density n in cm^-3? (Y/N) (default: N):" 
n
"Enter the conversion factor for number density n to cm^-3: (default: 1.0):" 

"Is the unit of temperature T in K? (Y/N) (default: N):" 
y
```
4. Set frequency range:
```julia
"Frequency range start (MHz) (default: 115):"

"Frequency range end (MHz) (default: 175):" 
150
"Frequency resolution (MHz) (default: 0.2):" 
0.183
```
5. Define Box size:
```julia
"Side of the Box size (pc), please give a Float (default: 50.0):" 
1000.0
"Side of the Box size (pixel) (default: 256):"
512
```
6. Enable Faraday rotation or not
```julia
"Do you want to include Faraday rotation in the computation of Q and U? (Y/N) (default: N):" 
y
"Faraday depth range start (rad/m^2) (default: -20):" 
-50
"Faraday depth range end (rad/m^2) (default: 20):" 
50
"Faraday depth resolution (rad/m^2) (default: 0.1):" 
0.5
```
7. Instrumental options
```julia
"Do you want to perform filtering for Synchrotron data? (Y/N) (default: N):" 
n
"Do you want to process all lines of sight (x, y, z), or choose specific ones? (Enter 'all' or 'choose') (default: All):"
all
```
8. Load emissivity file, see next section
```
"Enter the path to the interpolation file (default: /path/to/default/emissivity.dat):" 
/path/to/emissivity.dat
```
9. Apply **Wolfire et al. 2003** ([DOI:  
10.1086/368016](https://ui.adsabs.harvard.edu/abs/2003ApJ...587..278W/abstract)) electron density prescription
```julia
"Do you want to use the Wolfire et al. 2003 prescription? (Y/N) (default: N):" 
y
"Please enter the values for the constants:"
zeta (ionization rate by Cosmic Rays) (default: 1.8e-17): 
5e-16
"Geff (effective radiation field) (default: 1.0):" 

"omegaPAH (PAH grain alignment efficiency) (default: 0.5):" 

"XC (Conversion factor of H into C) (default: 0.00014):" 
```
### **Notes**
If no value is entered, the default value is automatically selected.
The session processes the selected simulations, using the defined parameters, and prepares the results for analysis.

### **Synchrotron Emissivity Interpolation**

The `EmissInterp.jl` module calculates synchrotron emissivities based on magnetic field strength (`BField`) and frequency (`nu`). It uses equations from **Padovani et al. 2021** ([DOI: 10.1051/0004-6361/202140799](https://doi.org/10.1051/0004-6361/202140799)) to compute the parallel and perpendicular components of emissivity through numerical integration.

The function `EmissInterp` takes arrays of magnetic field strengths (in microGauss) and frequencies (in MHz) as inputs. It calculates the emissivities for each combination of these values and writes the results to a file called `emissivity.dat`. The output file includes the magnetic field strength, frequency, parallel emissivity, and perpendicular emissivity in a tab-separated format.

This data can be directly used in **MOOSE** to analyze synchrotron emission for specific magnetic fields and frequencies of interest.

1. Define the ranges for magnetic field strengths and frequencies:
   ```julia
   BArray = [1.0, 5.0, 10.0]   # Magnetic fields in µG
   nuArray = [100, 200, 300]   # Frequencies in MHz
   EmissInterp(BArray, nuArray)
   ```
2. Use the generated emissivity.dat file in MOOSE for further processing and visualization.

---
## **Caveats**

While MOOSE provides powerful tools for processing synchrotron emission data, there are some important limitations and assumptions to keep in mind:

1. **Optically Thin Approximation**:
   - The calculations assume an optically thin medium, which is appropriate for high Galactic latitude regions. This means absorption effects are **not accounted for** in the current implementation.

2. **Simplistic Instrumental Model**:
   - The instrumental model used in MOOSE is very basic and does not capture the full complexity of real interferometers. Users should be cautious when interpreting results in scenarios that require detailed instrumental simulations.
---

## **Contributors**
- **Jack Berat** - Main Developer
