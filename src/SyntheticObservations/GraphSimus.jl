"""
    GraphStatsSimu(simu, list_LOS, dist_array, phi_array, vel_array, 
                   zeta, Geff, omegaPAH, XC, TCNM, TWNM, cmap)

Processes a given simulation directory to generate histograms and heatmaps for various data parameters.

# Arguments
- `simu`: Path to the simulation directory.
- `list_LOS`: List of lines of sight (e.g., ["x", "y", "z"]).
- `dist_array`: Array representing the distance values.
- `phi_array`: Array representing the Faraday depth values.
- `vel_array`: Array representing the velocity values.
- `zeta`: Parameter for Wolfire constants.
- `Geff`: Parameter for Wolfire constants.
- `omegaPAH`: Parameter for Wolfire constants.
- `XC`: Parameter for Wolfire constants.
- `TCNM`: Temperature of the cold neutral medium.
- `TWNM`: Temperature of the warm neutral medium.
- `cmap`: Colormap to be used for heatmaps.

# Returns
- None
"""
function GraphStatsSimu(simu, list_LOS, dist_array, phi_array, vel_array, zeta, Geff, omegaPAH, XC, TCNM, TWNM, cmap)

    println("------------------------------------------------------------------------------------------------")
    println("Processing simulation: $simu")
    
    # Read the simulation data
    println("Reading simulation data...")
    B1, B2, BLOS, V1, V2, VLOS, T, n = ReadSimulation(simu, "x", 1, 1, 1, 1e3)

    # Create histogram path
    path_histo = joinpath(simu, "Histograms")
    mkpath(path_histo)
    println("Output directory for histograms: $path_histo")
    
    # Calculate statistics for each parameter
    println("Calculating statistics for B1, B2, BLOS and Btot...")
    Btotal = Btot(B1, B2, BLOS)
    maxB1, indmaxB1, minB1, indminB1, meanB1, stdB1, skewB1, kurtB1 = CalculateStatistics(B1)
    maxB2, indmaxB2, minB2, indminB2, meanB2, stdB2, skewB2, kurtB2 = CalculateStatistics(B2)
    maxBLOS, indmaxBLOS, minBLOS, indminBLOS, meanBLOS, stdBLOS, skewBLOS, kurtBLOS = CalculateStatistics(BLOS)
    maxBtot, indmaxBtot, minBtot, indminBtot, meanBtot, stdBtot, skewBtot, kurtBtot = CalculateStatistics(Btotal)

    println("Calculating statistics for V1, V2, and VLOS...")
    Vtotal = Btot(V1, V2, VLOS)
    maxV1, indmaxV1, minV1, indminV1, meanV1, stdV1, skewV1, kurtV1 = CalculateStatistics(V1)
    maxV2, indmaxV2, minV2, indminV2, meanV2, stdV2, skewV2, kurtV2 = CalculateStatistics(V2)
    maxVLOS, indmaxVLOS, minVLOS, indminVLOS, meanVLOS, stdVLOS, skewVLOS, kurtVLOS = CalculateStatistics(VLOS)
    maxVtot, indmaxVtot, minVtot, indminVtot, meanVtot, stdVtot, skewVtot, kurtVtot = CalculateStatistics(Vtotal)

    println("Calculating statistics for temperature and density...")
    maxT, indmaxT, minT, indminT, meanT, stdT, skewT, kurtT = CalculateStatistics(T)
    maxn, indmaxn, minn, indminn, meann, stdn, skewn, kurtn = CalculateStatistics(n)
    
    # Generate histograms for magnetic field components
    println("Generating histogram for magnetic field components...")
    plot_histogram([B1, B2, BLOS], ["By", "Bz", "Bx"], L"B \, [\mu G]", "counts", 1000, joinpath(path_histo, "histogram_B.pdf"), (-15, 15), (0, 2e5))
    plot_histogram([Btotal], ["Btot"], L"B \, [\mu G]", "counts", 1000, joinpath(path_histo, "histogram_B_tot.pdf"), (-15, 15), (0, 2e5))
    # Generate histograms for velocity components
    println("Generating histogram for velocity components...")
    plot_histogram([V1, V2, VLOS], ["Vy", "Vz", "Vx"], L"V \, [km \, s^{-1}]", "counts", 1000, joinpath(path_histo, "histogram_V.pdf"), (-30, 30), (0, 1e5))
    plot_histogram([Vtotal], ["Vtot"], L"V \, [km \, s^{-1}]", "counts", 1000, joinpath(path_histo, "histogram_V_tot.pdf"), (-30, 30), (0, 2e5))
    # Generate histogram for log10 of density
    println("Generating histogram for log10 of density...")
    plot_histogram([log10.(n)], [L"\log_{10}(n)"], L"\log_{10}(n) \, [cm^{-3}]", "counts", 100, joinpath(path_histo, "histogram_log10n.pdf"))
    # Generate histogram for temperature
    println("Generating histogram for temperature...")
    plot_histogram([T], ["T"], L"T \, [K]", "counts", 100, joinpath(path_histo, "histogram_T.pdf"))

    # Process each line of sight
    for LOS in list_LOS
        path = joinpath(simu, LOS, "column_density")
        mkpath(path)
        println("Processing line of sight: $LOS")
        println("Output directory for column density: $path")
        
        # Read column density data from FITS files
        println("Reading column density data...")
        NHI = read_FITS_file(joinpath(simu, LOS, "HI", "NHI.fits"))
        NCNM = read_FITS_file(joinpath(simu, LOS, "HI", "NCNM.fits"))
        NLNM = read_FITS_file(joinpath(simu, LOS, "HI", "NLNM.fits"))
        NWNM = read_FITS_file(joinpath(simu, LOS, "HI", "NWNM.fits"))
        
        # Generate heatmaps for each column density
        println("Generating heatmaps for column density data...")
        plot_heatmap(NWNM, dist_array, dist_array, joinpath(path, "NWNM.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(NLNM, dist_array, dist_array, joinpath(path, "NLNM.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(NCNM, dist_array, dist_array, joinpath(path, "NCNM.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(NHI, dist_array, dist_array, joinpath(path, "NHI.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")

        # Calculate and plot gas fraction maps
        println("Calculating and generating gas fraction maps...")
        fCNMmassMap, fLNMmassMap, fWNMmassMap = GasFractionMap(n, T; TCNM=TCNM, TWNM=TWNM)
        plot_heatmap(fCNMmassMap, dist_array, dist_array, joinpath(path, "fCNM.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        fCNMVolumeMap, fLNMVolumeMap, fWNMVolumeMap = VolumeFractionMap(n, T; TCNM=TCNM, TWNM=TWNM)
        plot_heatmap(fCNMVolumeMap, dist_array, dist_array, joinpath(path, "fVolumeCNM.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")

        # Read rotation measure map
        println("Reading rotation measure map...")
        RMmap = read_FITS_file(joinpath(simu, LOS, "Synchrotron", "WithFaraday", "RMmap.fits"))
        plot_heatmap(RMmap, dist_array, dist_array, joinpath(path, "RM.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_histogram([RMmap], ["RM"], L"RM \, [rad m^{-2}]", "counts", 1000, joinpath(path_histo, "histogram_RM_$LOS.pdf"))

        # Calculate electron density using Wolfire's model
        println("Generating electron density histogram...")
        ne = read_FITS_file(joinpath(simu, LOS, "Synchrotron", "ne.fits"))
        plot_histogram([log10.(ne)], ["ne"], L"ne \, [cm^{-3}]", "counts", 1000, joinpath(path_histo, "histogram_ne_$LOS.pdf"))

        # Read FDF, HI, CNM, and WNM temperature-brightness maps
        println("Reading temperature-brightness maps for FDF, HI, CNM, and WNM...")
        FDF = read_FITS_file(joinpath(simu, LOS, "Synchrotron", "WithFaraday", "FDF.fits"))
        HI = read_FITS_file(joinpath(simu, LOS, "HI", "TbHI.fits"))
        CNM = read_FITS_file(joinpath(simu, LOS, "HI", "TbCNM.fits"))
        WNM = read_FITS_file(joinpath(simu, LOS, "HI", "TbWNM.fits"))

        # Calculate and plot moments for each map
        println("Calculating and generating moments for each map...")
        M0F, M1F, M2F = moments_map(FDF, phi_array)
        M0HI, M1HI, M2HI = moments_map(HI, vel_array)
        M0CNM, M1CNM, M2CNM = moments_map(CNM, vel_array)
        M0WNM, M1WNM, M2WNM = moments_map(WNM, vel_array)

        path_moments = joinpath(simu, LOS, "moments")
        mkpath(path_moments)
        println("Output directory for moments: $path_moments")
        
        plot_heatmap(M0F, dist_array, dist_array, joinpath(path_moments, "M0_Faraday.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(M1F, dist_array, dist_array, joinpath(path_moments, "M1_Faraday.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(M2F, dist_array, dist_array, joinpath(path_moments, "M2_Faraday.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(M0HI, dist_array, dist_array, joinpath(path_moments, "M0_HI.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(M1HI, dist_array, dist_array, joinpath(path_moments, "M1_HI.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(M2HI, dist_array, dist_array, joinpath(path_moments, "M2_HI.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(M0CNM, dist_array, dist_array, joinpath(path_moments, "M0_CNM.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(M1CNM, dist_array, dist_array, joinpath(path_moments, "M1_CNM.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(M2CNM, dist_array, dist_array, joinpath(path_moments, "M2_CNM.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(M0WNM, dist_array, dist_array, joinpath(path_moments, "M0_WNM.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(M1WNM, dist_array, dist_array, joinpath(path_moments, "M1_WNM.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")
        plot_heatmap(M2WNM, dist_array, dist_array, joinpath(path_moments, "M2_WNM.pdf"), colormap=cmap, xlabel="Distance (pc)", ylabel="Distance (pc)")

        # Generate histograms for the moments
        println("Generating histograms for moments...")
        plot_histogram([M0F], ["M0"], L"M0 \, [K rad m^{-2}]", "counts", 1000, joinpath(path_histo, "histogram_M0_Faraday_$LOS.pdf"))
        plot_histogram([M1F], ["M1"], L"M1 \, [rad m^{-2}]", "counts", 1000, joinpath(path_histo, "histogram_M1_Faraday_$LOS.pdf"))
        plot_histogram([M2F], ["M2"], L"M2 \, [rad m^{-2}]", "counts", 1000, joinpath(path_histo, "histogram_M2_Faraday_$LOS.pdf"))
        plot_histogram([M0HI], ["M0"], L"M0 \, [K km s^{-1}]", "counts", 1000, joinpath(path_histo, "histogram_M0_HI_$LOS.pdf"))
        plot_histogram([M1HI], ["M1"], L"M1 \, [km s^{-1}]", "counts", 1000, joinpath(path_histo, "histogram_M1_HI_$LOS.pdf"))
        plot_histogram([M2HI], ["M2"], L"M2 \, [km s^{-1}]", "counts", 1000, joinpath(path_histo, "histogram_M2_HI_$LOS.pdf"))
        plot_histogram([M0CNM], ["M0"], L"M0 \, [K km s^{-1}]", "counts", 1000, joinpath(path_histo, "histogram_M0_CNM_$LOS.pdf"))
        plot_histogram([M1CNM], ["M1"], L"M1 \, [km s^{-1}]", "counts", 1000, joinpath(path_histo, "histogram_M1_CNM_$LOS.pdf"))
        plot_histogram([M2CNM], ["M2"], L"M2 \, [km s^{-1}]", "counts", 1000, joinpath(path_histo, "histogram_M2_CNM_$LOS.pdf"))
        plot_histogram([M0WNM], ["M0"], L"M0 \, [K km s^{-1}]", "counts", 1000, joinpath(path_histo, "histogram_M0_WNM_$LOS.pdf"))
        plot_histogram([M1WNM], ["M1"], L"M1 \, [km s^{-1}]", "counts", 1000, joinpath(path_histo, "histogram_M1_WNM_$LOS.pdf"))
        plot_histogram([M2WNM], ["M2"], L"M2 \, [km s^{-1}]", "counts", 1000, joinpath(path_histo, "histogram_M2_WNM_$LOS.pdf"))
    end
end

# Ask user for inputs with defaults
TCNM = ask_user("Enter the temperature of the CNM", 200)
TWNM = ask_user("Enter the temperature of the WNM", 2000)

# Get instrumental parameters for HI
zeta, Geff, omegaPAH, XC = WolfireConstants()
PixelLength_pc, PixelLength_cm, BoxLength_pc, DistanceArray = DistanceParameters()
phiArray = FaradayParameters()
velArray = VelocityParameters()
cmap = :plasma

simu_list = split(read(`find /Users/jb270005/Desktop/simu_RAMSES -name "256"`, String), "\n")
pop!(simu_list)
sort!(simu_list)

list_LOS = ["x", "y", "z"]

for simu in simu_list
    GraphStatsSimu(simu, list_LOS, DistanceArray, phiArray, velArray, zeta, Geff, omegaPAH, XC, TCNM, TWNM, cmap)
end
