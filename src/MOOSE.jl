module MOOSE

using Crayons
using Dates
using JSON
using CSV
using DataFrames
using FITSIO
using StatsBase
using Interpolations
using Dierckx
using FFTW
using LinearAlgebra

include(joinpath("Utils", "ArrayMath.jl"))
include(joinpath("Utils", "Prompts.jl"))
include(joinpath("Utils", "Progress.jl"))
include(joinpath("Utils", "Errors.jl"))
include(joinpath("Utils", "InputValidation.jl"))

include(joinpath("FileIO", "FITSUtils.jl"))
include(joinpath("FileIO", "Header.jl"))
include(joinpath("FileIO", "ReadSimulation.jl"))
include(joinpath("FileIO", "WriteDataOnDisk.jl"))

include(joinpath("PhysicalParameters", "ConversionJyBeamtoK.jl"))
include(joinpath("PhysicalParameters", "BrightnessTemperature.jl"))
include(joinpath("PhysicalParameters", "Borientation.jl"))
include(joinpath("PhysicalParameters", "ElectronDensity.jl"))
include(joinpath("PhysicalParameters", "IntrinsicAngle.jl"))
include(joinpath("PhysicalParameters", "MagneticField.jl"))
include(joinpath("PhysicalParameters", "PolarizationAngle.jl"))
include(joinpath("PhysicalParameters", "PolarizationFraction.jl"))
include(joinpath("PhysicalParameters", "Pressure.jl"))
include(joinpath("PhysicalParameters", "RM.jl"))
include(joinpath("PhysicalParameters", "Constants.jl"))

include(joinpath("Frequencies", "FreqFile.jl"))

include(joinpath("Synchrotron", "EmissInterp.jl"))
include(joinpath("Synchrotron", "Pnu.jl"))
include(joinpath("Synchrotron", "ProcessSynchrotron.jl"))
include(joinpath("Synchrotron", "QUnu.jl"))
include(joinpath("Synchrotron", "Tnu.jl"))

include(joinpath("Faraday", "FaradayParameters.jl"))
include(joinpath("Faraday", "RMSynthesis.jl"))

include(joinpath("Filtering", "Filter.jl"))

include(joinpath("Statistics", "EffectiveWidth.jl"))
include(joinpath("Statistics", "Moments.jl"))
include(joinpath("Statistics", "RMS.jl"))
include(joinpath("Statistics", "PowerSpectrum.jl"))
include(joinpath("Statistics", "Statistics.jl"))

include(joinpath("SyntheticObservations", "SimulationDiscovery.jl"))
include(joinpath("SyntheticObservations", "DictHeaderParameters.jl"))
include(joinpath("SyntheticObservations", "InstrumentalParameters.jl"))
include(joinpath("SyntheticObservations", "MOOSE.jl"))
include(joinpath("SyntheticObservations", "MOOSE_from_config.jl"))

using .MOOSEFromConfig: MOOSE_from_config

export run_moose, MOOSE_from_config, MooseError, cli_error, config_error

end
