sigm(x, l, xh) = x > 0 ? (1/(1+exp(-l*(x-xh)))) : 0

H_inf(v, l, h) = v != 0 ? 1/(1+exp((v+l)/h)) : 0
#H_inf(v, l, h) = 1/(1+exp((v+l)/h))

B_inf(v, l, h) = v > 0 ? 1/(1+exp((v+l)/h)) : 0

Stim(t, stim_start, stim_end, photon_flux) = stim_start <= t <= stim_end ? photon_flux : 0.0

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

    #CONSTANTS
    RMAX = 1.0
    TMAX = 1.0
    PMAX = 1.0
    GMAX = 1.0

    #Open parameters
    (kP1, kP2, kP3, kP4, kP5, kP6,
    l1, h1, JMAX, HMAX,
    τR, τT, τP, τG, τJ, τH) = p

    @. dR = kP1*Stim(t, stim_start, stim_end, photon_flux)*(1-R/RMAX) - R/τR
    @. dT = kP2*R*(1-T/TMAX)- T/τT
    @. dP = kP3*T*(1-P/PMAX) - P/τP
    @. dG = -kP4*P*(1-G/-GMAX) - G/τG # Non-linear degradation
    @. dJ = ( kP5*G*(1-J/-JMAX) - J)/τJ
    @. dH = (kP6*H_inf(J, l1, h1)^2*(1-H/ HMAX) - H)/τH #to add this or not *H_inf(J, l1, h1)

    @. dA = (J+H) - A
    return nothing
end

function erg_ode!(du, u, p, t; stim_start = 0.0, stim_end = 1.0, photon_flux = 400.0)
    dB = view(du, 8)
    dM = view(du, 9)
    dC = view(du, 10)

    dO1 = view(du, 11)
    dO2 = view(du, 12)
    dO3 = view(du, 13)
    dO = view(du, 14)
    
    dERG = view(du, 15)
    
    A = view(u, 7)
    B = view(u, 8)
    M = view(u, 9)
    C = view(u, 10)

    O1 = view(u, 11)
    O2 = view(u, 12)
    O3 = view(u, 13)
    O = view(u, 14)
    
    ERG = view(u, 15)
    
    #Open parameters
    (photo_p..., #Splat the first several params away
    l2, h2,
    k1, k2, k3, k4, k5, k6, k7, 
    l3, h3, 
    τB, τM, τC, τO) = p
    
    #Run phototransduction first
    phototransduction_ode!(du, u, photo_p, t; stim_start = stim_start, stim_end = stim_end, photon_flux = photon_flux)

    @. dB = (-k1*B_inf(A, l2, h2) - B)/τB #This is transfer from PC to BPC #*(-A^3) This nonlinear term may be taking away
    @. dM = (k2*A - k3*B - M)/τM
    @. dC = (-k4*A - C)/τC 
    
    # Oscillator equations with delayed forcing applied to the O1 system:
    @. dO1 = O2
    @. dO2 = (-k5*O1 - k6*O2) + sigm(B, l3, h3)
    @. dO3 = (-k7*B*abs(O2) - O3)/τO
    @. dO = (O3 * O2) - O

    @. dERG = (A + B + M + C + O) - ERG
    nothing
end

function make_model(data, params; ex_vivo = true, stim_start = 0.0, stim_end = 1.0, photon_flux = 400.0)
    tspan = (0.0, data.t[end])
    u0 = [
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 
            0.0
        ]

    
    if ex_vivo
        params[22] = 0.0 #Silence k4 
        params[25] = 0.0 #Silence k7
    end

    prob = ODEProblem(
        (du, u, p, t) -> erg_ode!(du, u, p, t; stim_start = stim_start, stim_end = stim_end, photon_flux = photon_flux), 
        u0, tspan, params
    )
    return prob
end

function simulate_model(data, params; stim_start = 0.0, stim_end = 1.0, photon_flux = 400.0)
    prob = make_model(data, params; stim_start = stim_start, stim_end = stim_end,photon_flux = photon_flux)
    sol = solve(prob, Tsit5(), saveat=data.t, tstops=[stim_start, stim_end])
	ERG_t = map(t -> sol(t)[15], data.t)
    return sol, ERG_t
end

#--- Loss Function ---
function loss_static(data, params; channel = 1, kwargs...)
    # Define simulation time points from your data
    expERG = getchannel(data, channel).data_array[1,:,1]

    sol, ERG_t = simulate_model(data, params; kwargs...)
    n_length = length(ERG_t)
    # Compute the sum-of-squared error compared to experimental data
    return sum((ERG_t .- expERG).^2)#/n_length
end

function loss_static_abm(data_a, data_ab, data_abm, full_params; channel = 3, kwargs...)
    # Define simulation time points from your data
    expERG_A = getchannel(data_a, channel).data_array[1,:,1]
    expERG_AB = getchannel(data_ab, channel).data_array[1,:,1]
    expERG_ABM = getchannel(data_abm, channel).data_array[1,:,1]

    sol, ERG_t = simulate_model(data_a, full_params; kwargs...)
    #Compute the loss for each model
    a_wave = map(t -> sol(t)[7], data_a.t)
    loss_a = sum((a_wave .- expERG_A).^2)#/length(a_wave)
    #println("Loss A: $loss_a")

    ab_wave = map(t -> sol(t)[8], data_ab.t) .+ a_wave
    loss_ab = sum((ab_wave .- expERG_AB).^2)#/length(ab_wave)
    #println("Loss AB: $loss_ab")

    abm_wave = map(t -> sol(t)[15], data_abm.t)
    loss_abm = sum((abm_wave .- expERG_ABM).^2)#/length(abm_wave)
    #println("Loss ABM: $loss_abm")

    sum_loss = loss_a + loss_ab + loss_abm
    #println("Total Loss: $sum_loss")
    return sum_loss
end

#DO this for the phototransduction model only
function make_model_photo(data, params; ex_vivo = true, kwargs...)
    tspan = (0.0, data.t[end])
    u0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    prob = ODEProblem(
        (du, u, p, t) -> phototransduction_ode!(du, u, p, t; kwargs...), 
        u0, tspan, params
    )
    return prob
end

function simulate_model_photo(data, params; kwargs...)
    prob = make_model_photo(data, params; kwargs...)
    sol = solve(prob, Tsit5(), saveat=data.t, tstops=[stim_start, stim_end])
	ERG_t = map(t -> sol(t)[7], data.t)
    return sol, ERG_t
end

function loss_graded(data_series, params; channel = 1, stim_start = 0.0, stim_end = 1.0)
    # Define simulation time points from your data
    loss = 0.0
    ir_loss = 0.0
    weight_loss = 0.0

    loss_vals = []
    ir_loss_vals = []
    weight_loss_vals = []
    for (i, (k, data)) in enumerate(data_series)
        expERG = getchannel(data, channel).data_array[1,:,1]
        sol, ERG_t = simulate_model(data, params; stim_start, stim_end, photon_flux = k)
        a_wave = map(t -> sol(t)[7], data.t)
        
        #Compute the IR error
        loss_ir = abs(minimum(a_wave)) .- (minimum(expERG))
        ir_loss += loss_ir
        push!(loss_vals, loss_ir)
        
        # Compute the sum-of-squared error compared to experimental data
        n_length = length(ERG_t)
        loss_val = sum((a_wave .- expERG).^2)/n_length
        loss += loss_val
        push!(loss_vals, loss_val)

        #Compute the weighted loss
        weights = map(t -> exp(-t/200), data.t)  # Exponential decay 
        weight_loss_val = sum(weights .* (a_wave .- expERG).^2)/n_length
        weight_loss += weight_loss_val
        push!(weight_loss_vals, weight_loss_val)
    end
    return ir_loss, loss, weight_loss, ir_loss_vals, loss_vals, weight_loss_vals
    #return loss_vals
end