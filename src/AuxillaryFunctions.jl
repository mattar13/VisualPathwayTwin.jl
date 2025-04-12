#The stimulus function determines how much light is present at a given time.
Stim(t, stim_start, stim_end, photon_flux) = stim_start <= t <= stim_end ? photon_flux : 0.0

sigm(x, l, xh) = x > 0 ? (1/(1+exp(-l*(x-xh)))) : 0

H_inf(v, l, h) = v != 0 ? 1/(1+exp((v+l)/h)) : 0
#H_inf(v, l, h) = 1/(1+exp((v+l)/h))

B_inf(v, l, h) = v > 0 ? 1/(1+exp((v+l)/h)) : 0

#Below is a work in progress for a more complicated photoreceptor model
J(cGMP, Jmax) = (Jmax*cGMP^3)/(cGMP^3 + 10^3) #CNG currents
gL = 0.35 #Leakage conductance
EL = -77.0 #Leakage potential
ILeak(V) = gL*(V-EL)
Iphoto(V, cGMP, Jmax) = -J(cGMP, Jmax)*(1.0 - exp((V-8.5)/17.0))