function RobustSigma(Y; ZERO=false)
    # translation of the IDL function robust_sigma.pro
    # MAMD - 19/10/2023
    
    EPS = 1.0E-20
    if ZERO
        Y0 = 0.0
    else
        Y0 = median(Y)
    end

    MAD = median(abs.(Y .- Y0)) / 0.6745

    if MAD < EPS
        MAD = mean(abs.(Y .- Y0)) / 0.80
    end

    if MAD < EPS
        return 0.0
    end

    U = (Y .- Y0) ./ (6.0 * MAD)
    UU = U.^2
    Q = findall(UU .≤ 1.0)
    
    if length(Q) < 3
        println("ROBUST_SIGMA: This distribution is TOO WEIRD! Returning -1")
        return -1.0
    end

    N = count(!isnan, Y)
    NUMERATOR = sum((Y[Q] .- Y0).^2 .* (1.0 .- UU[Q]).^4)
    DEN1 = sum((1.0 .- UU[Q]) .* (1.0 .- 5.0 .* UU[Q]))
    SIGGMA = N * NUMERATOR / (DEN1 * (DEN1 - 1.0))

    if SIGGMA > 0.0
        return sqrt(SIGGMA)
    else
        return 0.0
    end
end
