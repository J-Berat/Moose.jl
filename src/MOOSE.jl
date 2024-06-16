module MOOSE

################################# USINGS  #################################

using FITSIO, DataFrames, CSV  # Data IO
using Interpolations, Dierckx   # Interpolation packages
using LinearAlgebra, SpecialFunctions # Maths
using QuadGK # Integration
using StatsBase, Statistics # Statistics
using LaTeXStrings, GLMakie, CairoMakie # Plotting packages
using StringEncodings # Operations on strings
using FFTW # Fourier
using ImageFiltering # Filtering images
using Random # Random generator, probability

################################# INCLUDES ################################# 

# Reading simulation
include("File_IO/ReadSimulation.jl") # Get all arrays B,n,T,V from a simulation

# Utils
include("Utils/Utils.jl") # useful functions

# Synthetic observations
include("SyntheticObservations/InstrumentalParameters.jl") # ask user for nu-,phhi and velArray
include("SyntheticObservations/DictHeaderParameters.jl") # Dictionnary for Header for synthetic obs
include("SyntheticObservations/Get3rdAxisValue.jl") # get Faraday Depth from the frame of Faraday cube
include("SyntheticObservations/MOOSE.jl") # Mock Observation Of Synchrotron Emission code

# Filtering
include("Filtering/Filter.jl") # Filter image by LOFAR filtering

# Statistics
include("Statistics/Statistics.jl") # compute basic statistics mean, median, std...
include("Statistics/Moments.jl") # compute moments from data

# Synchrotron
include("Synchrotron/EmissInterp.jl") # compute all relevant synchrotron equations
include("Synchrotron/Tnu.jl") # compute Tnu synchrotron
include("Synchrotron/QUnu.jl") # compute QUnu synchrotron
include("Synchrotron/Pnu.jl") # compute Pnu synchrotron
include("Synchrotron/ProcessSynchrotron.jl") # Process the whole synchrotron mock observation
include("Frequencies/FreqFile.jl")

# Physical Parameters
include("PhysicalParameters/Energies.jl") # compute Kinetic and Magnetic energies in CGS units
include("PhysicalParameters/RM.jl") # compute Rotation Measure in rad/m^2
include("PhysicalParameters/MagneticField.jl") # compute total and perpendicular to LOS magnetic field
include("PhysicalParameters/Pressure.jl") # compute Pressure from n and T
include("PhysicalParameters/ElectronDensity.jl") # compute electron density constant/prop_to_nH/Wolfire2003
include("PhysicalParameters/IntrinsicAngle.jl") # compute intrinsic polarization angle
include("PhysicalParameters/BrightnessTemperature.jl") # compute Brightness Temperature a-la single-dish
include("PhysicalParameters/PolarizationFraction.jl") # compute polarization fraction
include("PhysicalParameters/PolarizationAngle.jl") # compute polarization angle

# RM-Synthesis
include("Faraday/RMSynthesis.jl") # compute RM-Synthesis from Qnu and Unu
include("Faraday/FaradayWidth.jl") # compute the second order moment of the mean spectrum of a whole Faraday cube

# Writing files on disk
include("File_IO/Header.jl") # write header for HI, Faraday or Synchrotron cubes
include("File_IO/WriteDataOnDisk.jl") # Write HI, Faraday or Synchrotron cube on disk

################################# CONSTANTS ################################# 

# Fundamental Constants in IS
const ELECTRON_ENERGY_AT_REST_eV = 0.511e6 # eV electron energy at rest
const EMIN = 1e6 # minimum eV energy on the energy spectrum Padovani2021
const EMAX = 1e11 # maximum eV energy on the energy spectrum Padovani2021
const E_0 = 710e6 # eV Padovani 2021 SKA cf eq6
const C_m = 2.99792458e8 # speed of light in cm.s^-1

# in CGS
const E_CHARGE = 4.8032e-10 # electron charge in cm^3/2 g^1/2 s^-1 (esu-CGS)
const M_e = 9.109e-28 # electron mass in g
const M_p = 1.6726231e-24 # proton mass in g
const C = 2.99792458e10 # speed of light in cm.s^-1
const K_B = 1.380649e-16 # Boltzmann constant in cm^2 g s^-2 K^-1
const J_0 = 2.1e18 # prefactor of j_e en e^-1 s^-1 cm^-2 sr^-1, Padovani,Galli 2018 SKA
const RM_PREFACTOR = -0.81 #microG^-1 pc^-1 cm^-3

# Unit conversion
const EV_TO_ERG = 1.60218e-12 # eV->erg conversion factor
const JY_TO_MJY = 1e-6 # Jansky to MegaJansky conversion factor
const Wm2Hz_to_Jy = 1e26 # W.m^2.Hz to Jansky conversion factor
const MILI_TO_MICRO = 1e3 # mili -> mu conversion factor
const PARSEC_TO_CM = 3.086e18 # pc -> cm conversion factor

# prefactors
const PRE_NU_C = 3 * E_CHARGE / (4 * pi * M_e * C) # prefactor = 4.19 MHz Padovani,Galli 2018 SKA cf eq 3
const PRE_P = (sqrt(3) * E_CHARGE^3) / (2 * M_e * C^2) # prefactor eq.2 Padovani+ 2021

end
