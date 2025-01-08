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
- Compute Stokes parameters (\(I\), \(Q\), \(U\)).
- Perform Rotation Measure Synthesis (RM Synthesis) for Faraday rotation studies.
- Generate frequency files for simulated observations.
- Handle instrumental parameters such as frequency resolution and simulation box sizes.
- Support for CGS physical units.

---

## **Installation**
1. Clone this repository:
    ```bash
    git clone https://github.com/username/moose.git
    cd moose
    ```
2. Ensure Julia is installed (version 1.8+ recommended).
3. Install the required dependencies in Julia:
    ```julia
    using Pkg
    Pkg.add([
        "FITSIO",          # Data IO
        "DataFrames",
        "CSV",
        "Interpolations",  # Interpolation packages
        "Dierckx",
        "LinearAlgebra",   # Mathematics
        "SpecialFunctions",
        "QuadGK",          # Integration
        "StatsBase",       # Statistics
        "Statistics",
        "KernelDensity",
        "StringEncodings", # String operations
        "FFTW",            # Fourier
        "ImageFiltering",  # Image processing
        "Random",          # Random generator
        "Distributions"    # Probability distributions
    ])
    ```
---

## **Usage**
The primary tool in MOOSE is the **`MOOSE()`** function, which guides users interactively through the processing of simulated data. To launch MOOSE, use:
```julia
include("moose.jl")  # Load the main file
MOOSE()              # Start the interactive tool

---

## **Contributors**
- **JB** - Main Developer 
