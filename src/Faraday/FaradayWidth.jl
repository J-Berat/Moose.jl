function FaradayWidth(FDFcube, PhiArray)
    meanFDF = MeanSpectrum(FDFcube)
    m2 = moments(meanFDF, x=PhiArray)[3]
    return m2
end