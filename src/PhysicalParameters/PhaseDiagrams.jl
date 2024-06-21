function phase_diagram_pressure(n, P)
    histogram2d(vec(log10.(n)), vec(log10.(P)), xlabel=L"\log(n) [cm^{-3}]", ylabel=L"\log(P) [K.cm^{-3}]",colorbar_scale=:log10, grid=false, cbartitle="counts",title = L"P(n)")
end

function phase_diagram_temperature(n, T)
    histogram2d(vec(log10.(n)), vec(log10.(T)), xlabel=L"\log(n) [cm^{-3}]", ylabel=L"\log(T) [K]",colorbar_scale=:log10, grid = false, cbartitle="counts",title = L"T(n)")
end

function phase_diagram_ionizationfraction(Xe, T)
    histogram2d(vec(log10.(T)), vec(log10.(Xe)), xlabel=L"\log(T)[K]", ylabel=L"\log(Xe) ",colorbar_scale=:log10, grid = false, cbartitle="counts",title = L"Xe(T)")
end