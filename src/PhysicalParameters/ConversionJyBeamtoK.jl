"""
    ConversionJyK(I, frame, nu, theta_arcsec)

Convert intensity from Janskys to Kelvins.

# Arguments
- `I`: Intensity in Janskys.
- `frame`: Frame of reference.
- `nu`: Frequency in Hertz.
- `theta_arcsec`: Angular size in arcseconds.

# Returns
- An array representing the converted intensity in Kelvins.
"""
ConversionJyK(I,frame,nu,theta_arcsec) = 1.222*1e3 .* I ./ (nu^2*theta_arcsec^2)