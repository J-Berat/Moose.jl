"""
    header_params(; naxis, ctype1="", ctype2="", ctype3="", cunit1="", cunit2="", cunit3="", bunit="")

Utility helper to build the dictionaries describing FITS header metadata.
It keeps the key set consistent across entries in `DictHeader` while allowing
each dataset to specify its own axis labels, units, and dimensionality.
"""
function header_params(; naxis, ctype1="", ctype2="", ctype3="", cunit1="", cunit2="", cunit3="", bunit="")
    Dict(
        "naxis" => naxis,
        "ctype1" => ctype1,
        "ctype2" => ctype2,
        "ctype3" => ctype3,
        "cunit1" => cunit1,
        "cunit2" => cunit2,
        "cunit3" => cunit3,
        "bunit" => bunit,
    )
end

const DictHeader = Dict(
    "T_nu" => header_params(naxis=3, ctype3="FREQ", cunit3="Hz", bunit="K"),
    "Qnu" => header_params(naxis=3, ctype3="FREQ", cunit3="Hz", bunit="K"),
    "Unu" => header_params(naxis=3, ctype3="FREQ", cunit3="Hz", bunit="K"),
    "Pnu" => header_params(naxis=3, ctype3="FREQ", cunit3="Hz", bunit="K"),
    "TbHI" => header_params(naxis=3, ctype3="VEL", cunit1="deg", cunit2="deg", cunit3="km/s", bunit="K"),
    "tauHI" => header_params(naxis=3, ctype3="VEL", cunit3="km/s", bunit="K"),
    "TbthinHI" => header_params(naxis=3, ctype3="VEL", cunit3="km/s", bunit="K"),
    "TbCNM" => header_params(naxis=3, ctype3="VEL", cunit3="km/s", bunit="K"),
    "tauCNM" => header_params(naxis=3, ctype3="VEL", cunit3="km/s", bunit="K"),
    "TbthinCNM" => header_params(naxis=3, ctype3="VEL", cunit3="km/s", bunit="K"),
    "TbLNM" => header_params(naxis=3, ctype3="VEL", cunit3="km/s", bunit="K"),
    "tauLNM" => header_params(naxis=3, ctype3="VEL", cunit3="km/s", bunit="K"),
    "TbthinLNM" => header_params(naxis=3, ctype3="VEL", cunit3="km/s", bunit="K"),
    "TbWNM" => header_params(naxis=3, ctype3="VEL", cunit3="km/s", bunit="K"),
    "tauWNM" => header_params(naxis=3, ctype3="VEL", cunit3="km/s", bunit="K"),
    "TbthinWNM" => header_params(naxis=3, ctype3="VEL", cunit3="km/s", bunit="K"),
    "ne" => header_params(naxis=3, ctype3="VEL", cunit3="pc", bunit="cm^-3"),
    "NHI" => header_params(naxis=2, bunit="cm^-2"),
    "nCNM" => header_params(naxis=3, ctype3="VEL", cunit3="pc", bunit="cm^-3"),
    "nLNM" => header_params(naxis=3, ctype3="VEL", cunit3="pc", bunit="cm^-3"),
    "nWNM" => header_params(naxis=3, ctype3="VEL", cunit3="pc", bunit="cm^-3"),
    "NCNM" => header_params(naxis=2, bunit="cm^-2"),
    "NLNM" => header_params(naxis=2, bunit="cm^-2"),
    "NWNM" => header_params(naxis=2, bunit="cm^-2"),
    "intne" => header_params(naxis=2, bunit="cm^-2"),
    "sigmane" => header_params(naxis=2, bunit="cm^-3"),
    "intBLOS" => header_params(naxis=2, bunit="muG cm^2"),
    "sigmaBLOS" => header_params(naxis=2, bunit="muG"),
    "intBtotal" => header_params(naxis=2, bunit="muG cm^2"),
    "sigmaBtotal" => header_params(naxis=2, bunit="muG"),
    "intBperp" => header_params(naxis=2, bunit="muG/cm^2"),
    "sigmaT" => header_params(naxis=2, bunit="K"),
    "Pmax" => header_params(naxis=2, bunit="K"),
    "Pnumax" => header_params(naxis=2, bunit="rad.m^{-2}"),
    "RMmap" => header_params(naxis=2, bunit="rad.m^{-2}"),
    "FDF" => header_params(naxis=3, ctype3="FAR-DEPTH", cunit3="rad/m^2", bunit="K"),
    "realFDF" => header_params(naxis=3, ctype3="FAR-DEPTH", cunit3="rad/m^2", bunit="K"),
    "imagFDF" => header_params(naxis=3, ctype3="FAR-DEPTH", cunit3="rad/m^2", bunit="K"),
    "M0_Faraday" => header_params(naxis=2, bunit="K.rad/m^2"),
    "M1_Faraday" => header_params(naxis=2, bunit="rad/m^2"),
    "M2_Faraday" => header_params(naxis=2, bunit="rad/m^2"),
    "M0_HI" => header_params(naxis=2, bunit="K.km/s"),
    "M1_HI" => header_params(naxis=2, bunit="km/s"),
    "M2_HI" => header_params(naxis=2, bunit="km/s"),
)
#ListDataName = ["T_nu", "Qnu", "Unu", "TbHI", "tauHI", "TbthinHI", "TbCNM", "tauCNM", "TbthinCNM", "TbWNM", "tauWNM", "TbthinWNM", "ne",  "NHI", "NCNM", "NLNM", "NWNM", "intne", "intBLOS", "Pmax", "Pnumax", "RMmap", "FDF", "realFDF", "imagFDF"]

ListDataName = ["Qnu", "Unu", "ne", "intne", "intBLOS","Pnu", "Pmax", "Pnumax", "RMmap", "FDF", "realFDF", "imagFDF"]