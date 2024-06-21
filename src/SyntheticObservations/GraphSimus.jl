simu_list = split(read(`find /Users/jb270005/Desktop/simu_RAMSES -name "256"`,String),"\n")
pop!(simu_list)
simu_list = sort!(simu_list)

list_LOS = ["x", "y", "z"]

# Box size
BoxLength = 50 # pc
BoxResolution = 256 # pixel
PixelLength_pc = BoxLength/BoxResolution
PixelLength_cm = PixelLength_pc*PARSEC_TO_CM
Dstart = 0
Dend = 50
dD = PixelLength_pc
DistanceArray = range(start=Dstart,stop=Dend,step=dD)

TCNM = 200
TWNM = 2000

# Faraday depth range in rad/m^2
Phistart = -10
Phiend = 10
dPhi = 0.25 # Faraday depth resolution
# velocity range in km/s
velstart = -30
velend = 30
dvel = 1.288214969124 # velocity resolution
# Distance range in pc 
Dstart = 0
Dend = 50
dD = PixelLength_pc

DistanceArray = range(start=Dstart,stop=Dend,step=dD) # pc
phiArray = range(start=Phistart,stop=Phiend,step=dPhi) # rad/m^2
velArray = range(start=velstart,stop=velend,step=dvel); # km/s

zeta, Geff, omegaPAH, XC = WolfireConstants()

for simu in simu_list[1:end]
    println("------------------------------------------------------------------------------------------------")
    println("The used simulation is the following : $simu")
    
    B1,B2,BLOS,V1,V2,VLOS,T,n = ReadSimulation(simu,"x", 1, 1, 1, 1e3)

    path_histo = joinpath(simu,"Histograms")
    mkpath(path_histo)
    
    maxB1, indmaxB1, minB1, indminB1, meanB1, stdB1, skewB1, kurtB1 = CalculateStatistics(B1)
    maxB2, indmaxB2, minB2, indminB2, meanB2, stdB2, skewB2, kurtB2 = CalculateStatistics(B2)
    maxBLOS, indmaxBLOS, minBLOS, indminBLOS, meanBLOS, stdBLOS, skewBLOS, kurtBLOS = CalculateStatistics(BLOS)
    maxV1, indmaxV1, minV1, indminV1, meanV1, stdV1, skewV1, kurtV1 = CalculateStatistics(V1)
    maxV2, indmaxV2, minV2, indminV2, meanV2, stdV2, skewV2, kurtV2 = CalculateStatistics(V2)
    maxVLOS, indmaxVLOS, minVLOS, indminVLOS, meanVLOS, stdVLOS, skewVLOS, kurtVLOS = CalculateStatistics(VLOS)
    maxT, indmaxT, minT, indminT, meanT, stdT, skewT, kurtT = CalculateStatistics(T)
    maxn, indmaxn, minn, indminn, meann, stdn, skewn, kurtn = CalculateStatistics(n)

    f = Figure(size = (1200, 800))
    ax = Axis(f[1,1], xlabel = L"B \, [\mu G]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
    hist!(ax, vec(B1), color=:red,bins=1000,label="By")
    hist!(ax, vec(B2), color=:blue,bins=1000,label="Bz")
    hist!(ax, vec(BLOS), color=:green,bins=1000,label="Bx")
    vlines!(ax,[meanB1],label = L"\langle B_z \rangle = \langle B_y \rangle",linestyle=:dash,color=:darkred)
    vlines!(ax,[meanBLOS], label = L"\langle B_x \rangle",linestyle=:dash,color=:darkgreen)
    xlims!(-15,15)
    ylims!(0,2e5)
    f[1,2] = Legend(f, ax, framevisible=false)
    save(joinpath(path_histo,"histogram_B.pdf"),f)

    f2 = Figure(size = (1200, 800))
    ax2 = Axis(f2[1,1], xlabel = L"V \, [km \, s^{-1}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
    hist!(ax2, vec(V1), color=:red,bins=1000,label="Vy")
    hist!(ax2, vec(V2), color=:blue,bins=1000,label="Vz")
    hist!(ax2, vec(VLOS), color=:green,bins=1000,label="Vx")
    vlines!(ax2,[meanV1],label = L"\langle V_z \rangle",linestyle=:dash,color=:darkred)
    vlines!(ax2,[meanV2],label = L"\langle V_y \rangle",linestyle=:dash,color=:darkblue)
    vlines!(ax2,[meanVLOS], label = L"\langle V_x \rangle",linestyle=:dash,color=:darkgreen)
    xlims!(-30,30)
    ylims!(0,1e5)
    f2[1,2] = Legend(f2, ax2, framevisible=false)
    save(joinpath(path_histo,"histogram_V.pdf"),f2)
	
    f3 = Figure(size = (1200, 800))
    ax3 = Axis(f3[1,1], xlabel = L"\log_{10}(n) \, [cm^{-3}]", ylabel = "counts", xgridvisible = false, ygridvisible = false, yscale=log10)
    hist!(ax3, vec(log10.(n)), color=:blue, bins=100,label=L"\log_{10}(n)")
    vlines!(ax3,[log10.(meann)],label = L"\langle \log_{10}(n) \rangle",linestyle=:dash,color=:black)
    f3[1,2] = Legend(f3, ax3, framevisible=false)
    save(joinpath(path_histo,"histogram_log10n.pdf"),f3)

    f4 = Figure(size = (1200, 800))
    ax4 = Axis(f4[1,1], xlabel = L"T \, [K]", ylabel = "counts", xgridvisible = false, ygridvisible = false, yscale=log10)
    hist!(ax4, vec(T), color=:blue, bins=100, weights = n)
    vlines!(ax4, [meanT],label = L"\langle T \rangle",linestyle=:dash,color=:black)
    save(joinpath(path_histo,"histogram_T.pdf"),f4)
    
    #  LOS
    for LOS in list_LOS
        path = joinpath(simu,LOS,"column_density")
    	mkpath(path)
        println("LOS: $LOS")
        NHI = read(FITS(joinpath(simu,LOS,"HI","NHI.fits"))[1])
        NCNM = read(FITS(joinpath(simu,LOS,"HI","NCNM.fits"))[1])
        NLNM = read(FITS(joinpath(simu,LOS,"HI","NLNM.fits"))[1])
        NWNM = read(FITS(joinpath(simu,LOS,"HI","NWNM.fits"))[1])
    	
    	fig = Figure()
    	ax = Axis(fig[1,1])
    	hm = heatmap!(ax, DistanceArray, DistanceArray, NWNM)
    	Colorbar(fig[:, end+1], hm)
    	save(joinpath(path,"NWNM.pdf"),fig)
    
    	fig2 = Figure()
    	ax2 = Axis(fig2[1,1])
    	hm = heatmap!(ax2, DistanceArray, DistanceArray, NLNM)
    	Colorbar(fig2[:, end+1], hm)
    	save(joinpath(path,"NLNM.pdf"),fig2)
    
    	fig3 = Figure()
    	ax3 = Axis(fig3[1,1])
    	hm = heatmap!(ax3, DistanceArray, DistanceArray, NCNM)
    	Colorbar(fig3[:, end+1], hm)
    	save(joinpath(path,"NCNM.pdf"),fig3)
    
    	fig4 = Figure()
    	ax4 = Axis(fig4[1,1])
    	hm = heatmap!(ax4, DistanceArray, DistanceArray, NHI)
    	Colorbar(fig4[:, end+1], hm)
    	save(joinpath(path,"NHI.pdf"),fig4)

        fCNMmassMap, fLNMmassMap, fWNMmassMap = GasFractionMap(n, T; TCNM=200, TWNM=2000)
    	fig15 = Figure()
    	ax15 = Axis(fig15[1,1])
    	hm = heatmap!(ax15, DistanceArray, DistanceArray, fCNMmassMap)
    	Colorbar(fig15[:, end+1], hm)
    	save(joinpath(path,"fCNM.pdf"),fig15)

        fCNMVolumeMap, fLNMVolumeMap, fWNMVolumeMap = VolumeFractionMap(n, T; TCNM=200, TWNM=2000)
    	fig16 = Figure()
    	ax16 = Axis(fig16[1,1])
    	hm = heatmap!(ax16, DistanceArray, DistanceArray, fCNMVolumeMap)
    	Colorbar(fig16[:, end+1], hm)
    	save(joinpath(path,"fVolumeCNM.pdf"),fig16)

        RMmap = read(FITS(joinpath(simu,LOS,"Synchrotron","WithFaraday","RMmap.fits"))[1])
    	fig16 = Figure()
    	ax16 = Axis(fig16[1,1])
    	hm = heatmap!(ax16, DistanceArray, DistanceArray, RMmap)
    	Colorbar(fig16[:, end+1], hm)
    	save(joinpath(path,"RM.pdf"),fig16)

        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"RM \, [rad m^{-2}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, vec(RMmap), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_RM_$LOS.pdf"),f)
        
        ne = Wolfire_ne(zeta, Geff, omegaPAH, XC, T, n) 
        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"ne \, [cm^{-3}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, vec(log10.(ne)), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_ne_$LOS.pdf"),f)
        
        FDF = read(FITS(joinpath(simu,LOS,"Synchrotron","WithFaraday","FDF.fits"))[1])
        HI = read(FITS(joinpath(simu,LOS,"HI","TbHI.fits"))[1])
        CNM = read(FITS(joinpath(simu,LOS,"HI","TbCNM.fits"))[1])
        WNM = read(FITS(joinpath(simu,LOS,"HI","TbWNM.fits"))[1])

        M0F = zeros(size(FDF,1),size(FDF,2))
        M1F = zeros(size(FDF,1),size(FDF,2))
        M2F = zeros(size(FDF,1),size(FDF,2))
        M0HI = zeros(size(HI,1),size(HI,2))
        M1HI = zeros(size(HI,1),size(HI,2))
        M2HI = zeros(size(HI,1),size(HI,2))
        M0CNM = zeros(size(CNM,1),size(CNM,2))
        M1CNM = zeros(size(CNM,1),size(CNM,2))
        M2CNM = zeros(size(CNM,1),size(CNM,2))
        M0WNM = zeros(size(WNM,1),size(WNM,2))
        M1WNM = zeros(size(WNM,1),size(WNM,2))
        M2WNM = zeros(size(WNM,1),size(WNM,2))
        
        CF = MapComplexity(FDF, phiArray)
        CHI = MapComplexity(HI, velArray)
        CCNM = MapComplexity(CNM, velArray)
        CWNM = MapComplexity(WNM, velArray)

        path_complexity = joinpath(simu,LOS,"Complexity")
    	mkpath(path_complexity)
        
        fig4 = Figure()
    	ax4 = Axis(fig4[1,1])
    	hm = heatmap!(ax4, DistanceArray, DistanceArray, CF)
    	Colorbar(fig4[:, end+1], hm)
    	save(joinpath(path_complexity, "complexity_Faraday.pdf"),fig4)

    	fig5 = Figure()
    	ax5 = Axis(fig5[1,1])
    	hm = heatmap!(ax5, DistanceArray, DistanceArray, CHI)
    	Colorbar(fig5[:, end+1], hm)
    	save(joinpath(path_complexity, "complexity_HI.pdf"),fig5)

    	fig5 = Figure()
    	ax5 = Axis(fig5[1,1])
    	hm = heatmap!(ax5, DistanceArray, DistanceArray, CCNM)
    	Colorbar(fig5[:, end+1], hm)
    	save(joinpath(path_complexity, "complexity_CNM.pdf"), fig5)

        fig5 = Figure()
    	ax5 = Axis(fig5[1,1])
    	hm = heatmap!(ax5, DistanceArray, DistanceArray, CWNM)
    	Colorbar(fig5[:, end+1], hm)
    	save(joinpath(path_complexity, "complexity_WNM.pdf"), fig5)
        
        
        for i in 1:size(n,1)
            for j in 1:size(n,2)
                m0, m1, m2 = moments(FDF[i,j,:]; x=phiArray)
                M0F[i,j] = m0
                M1F[i,j] = m1
                M2F[i,j] = m2
                m0, m1, m2 = moments(HI[i,j,:]; x=velArray)
                M0HI[i,j] = m0
                M1HI[i,j] = m1
                M2HI[i,j] = m2
                m0, m1, m2 = moments(CNM[i,j,:]; x=velArray)
                M0CNM[i,j] = m0
                M1CNM[i,j] = m1
                M2CNM[i,j] = m2
                m0, m1, m2 = moments(WNM[i,j,:]; x=velArray)
                M0WNM[i,j] = m0
                M1WNM[i,j] = m1
                M2WNM[i,j] = m2
            end
        end
        
        path_moments = joinpath(simu,LOS,"moments")
    	mkpath(path_moments)
    	fig4 = Figure()
    	ax4 = Axis(fig4[1,1])
    	hm = heatmap!(ax4, DistanceArray, DistanceArray, M0F)
    	Colorbar(fig4[:, end+1], hm)
    	save(joinpath(path_moments,"M0_Faraday.pdf"),fig4)

    	fig5 = Figure()
    	ax5 = Axis(fig5[1,1])
    	hm = heatmap!(ax5, DistanceArray, DistanceArray, M1F)
    	Colorbar(fig5[:, end+1], hm)
    	save(joinpath(path_moments,"M1_Faraday.pdf"),fig5)

    	fig5 = Figure()
    	ax5 = Axis(fig5[1,1])
    	hm = heatmap!(ax5, DistanceArray, DistanceArray, M2F)
    	Colorbar(fig5[:, end+1], hm)
    	save(joinpath(path_moments,"M2_Faraday.pdf"), fig5)

    	fig6 = Figure()
    	ax6 = Axis(fig6[1,1])
    	hm = heatmap!(ax6, DistanceArray, DistanceArray, M0HI)
    	Colorbar(fig6[:, end+1], hm)
    	save(joinpath(path_moments,"M0_HI.pdf"), fig6)

    	fig7 = Figure()
    	ax7 = Axis(fig7[1,1])
    	hm = heatmap!(ax7, DistanceArray, DistanceArray, M1HI)
    	Colorbar(fig7[:, end+1], hm)
    	save(joinpath(path_moments,"M1_HI.pdf"), fig7)

    	fig8 = Figure()
    	ax8 = Axis(fig8[1,1])
    	hm = heatmap!(ax8, DistanceArray, DistanceArray, M2HI)
    	Colorbar(fig8[:, end+1], hm)
    	save(joinpath(path_moments,"M2_HI.pdf"), fig8)

    	fig9 = Figure()
    	ax9 = Axis(fig9[1,1])
    	hm = heatmap!(ax9, DistanceArray, DistanceArray, M0CNM)
    	Colorbar(fig9[:, end+1], hm)
    	save(joinpath(path_moments,"M0_CNM.pdf"), fig9)

    	fig10 = Figure()
    	ax10 = Axis(fig10[1,1])
    	hm = heatmap!(ax10, DistanceArray, DistanceArray, M1CNM)
    	Colorbar(fig10[:, end+1], hm)
    	save(joinpath(path_moments,"M1_CNM.pdf"), fig10)

    	fig11 = Figure()
    	ax11 = Axis(fig11[1,1])
    	hm = heatmap!(ax11, DistanceArray, DistanceArray, M2CNM)
    	Colorbar(fig11[:, end+1], hm)
    	save(joinpath(path_moments,"M2_CNM.pdf"), fig11)
        
    	fig12 = Figure()
    	ax12 = Axis(fig12[1,1])
    	hm = heatmap!(ax12, DistanceArray, DistanceArray, M0WNM)
    	Colorbar(fig12[:, end+1], hm)
    	save(joinpath(path_moments,"M0_WNM.pdf"), fig12)
        
    	fig13 = Figure()
    	ax13 = Axis(fig13[1,1])
    	hm = heatmap!(ax13, DistanceArray, DistanceArray, M1WNM)
    	Colorbar(fig13[:, end+1], hm)
    	save(joinpath(path_moments,"M1_WNM.pdf"), fig13)

    	fig14 = Figure()
    	ax14 = Axis(fig14[1,1])
    	hm = heatmap!(ax14, DistanceArray, DistanceArray, M2WNM)
    	Colorbar(fig14[:, end+1], hm)
    	save(joinpath(path_moments,"M2_WNM.pdf"), fig14)

        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"M0 \, [K rad m^{-2}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, filter(!isnan, vec(M0F)), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_M0_Faraday_$LOS.pdf"),f)

        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"M1 \, [rad m^{-2}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, filter(!isnan, vec(M1F)), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_M1_Faraday_$LOS.pdf"),f)

        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"M2 \, [rad m^{-2}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, filter(!isnan, vec(M2F)), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_M2_Faraday_$LOS.pdf"),f)

        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"M0 \, [K km s^{-1}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, filter(!isnan, vec(M0HI)), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_M0_HI_$LOS.pdf"),f)

        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"M1 \, [km s^{-1}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, filter(!isnan, vec(M1HI)), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_M1_HI_$LOS.pdf"),f)

        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"M2 \, [km s^{-1}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, filter(!isnan, vec(M2HI)), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_M2_HI_$LOS.pdf"),f)

        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"M0 \, [K km s^{-1}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, filter(!isnan, vec(M0CNM)), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_M0_CNM_$LOS.pdf"),f)

        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"M1 \, [km s^{-1}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, filter(!isnan, vec(M1CNM)), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_M1_CNM_$LOS.pdf"),f)

        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"M2 \, [km s^{-1}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, filter(!isnan, vec(M2CNM)), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_M2_CNM_$LOS.pdf"),f)

        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"M0 \, [K km s^{-1}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, filter(!isnan, vec(M0WNM)), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_M0_WNM_$LOS.pdf"),f)

        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"M1 \, [km s^{-1}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, filter(!isnan, vec(M1WNM)), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_M1_WNM_$LOS.pdf"),f)

        f = Figure(size = (1200, 800))
        ax = Axis(f[1,1], xlabel = L"M2 \, [km s^{-1}]", ylabel = "counts", xgridvisible = false, ygridvisible = false)
        hist!(ax, vec(filter(!isnan, M2WNM)), color=:blue,bins=1000)
        save(joinpath(path_histo,"histogram_M2_WNM_$LOS.pdf"),f)
    end
end