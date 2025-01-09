# **MOOSE**
**Mock Observation Of Synchrotron Emission**

MOOSE is an interactive tool designed for processing simulated synchrotron emission data. It provides functionalities for calculating intensity, polarization, and Faraday rotation maps, while supporting various physical and instrumental parameters.

---

## **Table of Contents**
1. [Main Features](#main-features)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Contributors](#contributors)

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
9. Apply Wolfire et al. 2003 electron density prescription
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

## **Contributors**
- **Jack Berat** - Main Developer
