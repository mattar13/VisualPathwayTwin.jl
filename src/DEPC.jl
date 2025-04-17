function phototransduction_ode!(du, u, p, t; stim_start = 0.0, stim_end = 1.0, photon_flux = 400.0)
    #Extract the parameters
    dR = view(du, 1)
    dT = view(du, 2)
    dP = view(du, 3)
    dG = view(du, 4)
    dJ = view(du, 5)
    dH = view(du, 6)

    dA = view(du, 7)

    R = view(u, 1)
    T = view(u, 2)
    P = view(u, 3)
    G = view(u, 4)
    J = view(u, 5)
    H = view(u, 6)

    A = view(u, 7)

    #Open parameters
    (aC, kR1, kF2, kR2, kF3, kR3, kHYDRO, kREC, 
    G0, iDARK, kg, 
    C_m, V0, gREST, l1, h1, gH, τH) = p

    #CONSTANTS
    G0 = 4.0
    kg = 20.0
    iDARK = 5040.0

    Φ=Stim(t, stim_start, stim_end, photon_flux)
    @. dR = aC*Φ - kR1*R
    @. dT = kF2*R*(1-T) - kR2*T
    @. dP = kF3*T*(1-P) - kR3*P
    @. dG = -kHYDRO*P*G + kREC*(G0 - G) # Non-linear degradation

    @. dJ = -iDARK * J∞(G, kg) - J#( kP5*G*(1-J/-JMAX) - J)/τJ
    @. dH = (gH*H_inf(A, l1, h1)*(A-0.25) - H)/τH #0.0#(kP6*H_inf(J, l1, h1)^2*(1-H/ HMAX) - H)/τH #to add this or not *H_inf(J, l1, h1)
    
    R_m = 10
    @. dA = -((J+V0)/R_m + H + gREST*(A-0.0))/C_m 
    return nothing
end
