# ============================================================
# photoreceptor.jl - Rod and cone photoreceptor dynamics
# ============================================================

# ── 1. Default Parameters ───────────────────────────────────

"""
    default_rod_params()

Return default parameters for the rod photoreceptor model as a named tuple.
All rates are in 1/s and converted to 1/ms inside the model function.

State vector (19 vars): [V, mKv, hKv, mCa, mKCa, Ca_s, Ca_f, Cab_ls, Cab_hs,
                         Cab_lf, Cab_hf, Rh, Rhi, Tr, PDE, Ca_photo, Cab_photo, cGMP, Glu]
"""
function default_rod_params()
    return (
        # Membrane properties
        C_m = 0.02,              # nF — membrane capacitance
        g_L = 0.35,              # nS — leak conductance
        E_L = -77.0,             # mV — leak reversal potential

        # Light capture
        eta = 0.67,              # quantum efficiency (dimensionless)

        # Phototransduction cascade
        alpha1 = 50.0,           # 1/s — Rh activation rate
        alpha2 = 0.0003,         # 1/s — Rh* to Rhi reverse rate
        alpha3 = 0.03,           # 1/s — Rhi inactivation rate
        epsilon = 0.5,           # 1/(s*µM) — transducin activation by Rh*
        beta1 = 2.5,             # 1/s — transducin inactivation
        tau1 = 0.2,              # 1/(s*µM) — PDE activation by Tr*
        tau2 = 5.0,              # 1/s — PDE inactivation
        T_tot = 1000.0,          # µM — total transducin
        PDE_tot = 100.0,         # µM — total PDE

        # cGMP dynamics
        J_max = 5040.0,          # pA — max photocurrent
        b = 0.25,                # µM/(s*pA) — Ca influx per photocurrent
        gamma_Ca = 50.0,         # 1/s — Ca extrusion from outer segment
        C0 = 0.1,                # µM — basal Ca
        k1 = 0.2,                # 1/(s*µM) — Ca buffer on-rate (outer segment)
        k2 = 0.8,                # 1/s — Ca buffer off-rate (outer segment)
        eT = 500.0,              # µM — total Ca buffer (outer segment)
        A_max = 65.6,            # µM/s — max cyclase rate
        K_c = 0.1,               # µM — cyclase Ca half-saturation
        nu = 0.4,                # 1/s — basal cGMP hydrolysis
        sigma = 1.0,             # 1/(s*µM) — PDE-mediated cGMP hydrolysis

        # I_H (hyperpolarization-activated current)
        g_H = 1.5,               # nS — max I_H conductance
        E_H = -32.0,             # mV — I_H reversal potential
        V_h_half = -70.0,        # mV — I_H half-activation voltage
        k_h = -7.0,              # mV — I_H activation slope

        # Voltage-gated K+ (delayed rectifier)
        g_Kv = 2.0,              # nS — max Kv conductance
        E_K = -74.0,             # mV — K+ reversal potential

        # Voltage-gated Ca2+
        g_Ca = 0.7,              # nS — max Ca conductance
        Ca_o = 1600.0,           # µM — extracellular Ca

        # Ca-activated Cl-
        g_Cl = 2.0,              # nS — max Cl conductance
        E_Cl = -20.0,            # mV — Cl reversal potential

        # Ca-activated K+
        g_KCa = 5.0,             # nS — max KCa conductance

        # Ca exchangers at inner segment membrane
        J_ex_max = 20.0,         # pA — max exchanger current
        K_ex = 0.2,              # µM — exchanger half-saturation
        J_ex2_max = 20.0,         # pA — max secondary exchanger current
        K_ex2 = 0.5,             # µM — secondary exchanger half-saturation

        # Intracellular Ca dynamics (inner segment)
        F = 96485.33212,         # C/mol — Faraday constant
        V1 = 3.812e-13,          # L — submembrane volume
        D_Ca = 6.0e-8,           # cm²/s — Ca diffusion coefficient
        S1 = 3.142e-8,           # cm² — effective surface area
        V2 = 5.236e-13,          # L — bulk cytoplasmic volume
        delta = 3e-5,            # dimensionless — geometry factor

        # Ca buffers (inner segment)
        Lb1 = 2.0,               # 1/(s*µM) — low-affinity buffer on-rate
        Lb2 = 1.0,               # 1/s — low-affinity buffer off-rate
        Hb1 = 1.11,              # 1/(s*µM) — high-affinity buffer on-rate
        Hb2 = 1.0,               # 1/s — high-affinity buffer off-rate
        B_L = 500.0,             # µM — total low-affinity buffer
        B_H = 300.0,             # µM — total high-affinity buffer

        # Glutamate release
        alpha_Glu = 1.0,         # dimensionless — max glutamate release
        V_Glu_half = -40.0,      # mV — release half-activation voltage
        V_Glu_slope = 5.0,       # mV — release slope
        tau_Glu = 5.0            # ms — glutamate clearance time constant
    )
end

"""
    default_cone_params()

Return default parameters for the cone photoreceptor model as a named tuple.
This is a simplified 6-state model. FOR NOW THIS IS THE SAME AS RODS BUT WE CAN ADJUST

State vector (6 vars): [R*, G, Ca, V, h, Glu]
"""
function default_cone_params()
    return (
        # Membrane properties
        C_m = 0.02,              # nF — membrane capacitance
        g_L = 0.35,              # nS — leak conductance
        E_L = -77.0,             # mV — leak reversal potential

        # Light capture
        eta = 0.67,              # quantum efficiency (dimensionless)

        # Phototransduction cascade
        alpha1 = 50.0,           # 1/s — Rh activation rate
        alpha2 = 0.0003,         # 1/s — Rh* to Rhi reverse rate
        alpha3 = 0.03,           # 1/s — Rhi inactivation rate
        epsilon = 0.5,           # 1/(s*µM) — transducin activation by Rh*
        beta1 = 2.5,             # 1/s — transducin inactivation
        tau1 = 0.2,              # 1/(s*µM) — PDE activation by Tr*
        tau2 = 5.0,              # 1/s — PDE inactivation
        T_tot = 1000.0,          # µM — total transducin
        PDE_tot = 100.0,         # µM — total PDE

        # cGMP dynamics
        J_max = 5040.0,          # pA — max photocurrent
        b = 0.25,                # µM/(s*pA) — Ca influx per photocurrent
        gamma_Ca = 50.0,         # 1/s — Ca extrusion from outer segment
        C0 = 0.1,                # µM — basal Ca
        k1 = 0.2,                # 1/(s*µM) — Ca buffer on-rate (outer segment)
        k2 = 0.8,                # 1/s — Ca buffer off-rate (outer segment)
        eT = 500.0,              # µM — total Ca buffer (outer segment)
        A_max = 65.6,            # µM/s — max cyclase rate
        K_c = 0.1,               # µM — cyclase Ca half-saturation
        nu = 0.4,                # 1/s — basal cGMP hydrolysis
        sigma = 1.0,             # 1/(s*µM) — PDE-mediated cGMP hydrolysis

        # I_H (hyperpolarization-activated current)
        g_H = 1.5,               # nS — max I_H conductance
        E_H = -32.0,             # mV — I_H reversal potential
        V_h_half = -70.0,        # mV — I_H half-activation voltage
        k_h = -7.0,              # mV — I_H activation slope

        # Voltage-gated K+ (delayed rectifier)
        g_Kv = 2.0,              # nS — max Kv conductance
        E_K = -74.0,             # mV — K+ reversal potential

        # Voltage-gated Ca2+
        g_Ca = 0.7,              # nS — max Ca conductance
        Ca_o = 1600.0,           # µM — extracellular Ca

        # Ca-activated Cl-
        g_Cl = 2.0,              # nS — max Cl conductance
        E_Cl = -20.0,            # mV — Cl reversal potential

        # Ca-activated K+
        g_KCa = 5.0,             # nS — max KCa conductance

        # Ca exchangers at inner segment membrane
        J_ex_max = 20.0,         # pA — max exchanger current
        K_ex = 0.2,              # µM — exchanger half-saturation
        J_ex2_max = 5.0,         # pA — max secondary exchanger current
        K_ex2 = 0.5,             # µM — secondary exchanger half-saturation

        # Intracellular Ca dynamics (inner segment)
        F = 96485.33212,         # C/mol — Faraday constant
        V1 = 3.812e-13,          # L — submembrane volume
        V2 = 5.236e-13,          # L — bulk cytoplasmic volume
        S1 = 3.142e-8,           # cm² — effective surface area
        delta = 0.05,            # dimensionless — geometry factor
        D_Ca = 6.0e-8,           # cm²/s — Ca diffusion coefficient

        # Ca buffers (inner segment)
        Lb1 = 2.0,               # 1/(s*µM) — low-affinity buffer on-rate
        Lb2 = 1.0,               # 1/s — low-affinity buffer off-rate
        Hb1 = 1.11,              # 1/(s*µM) — high-affinity buffer on-rate
        Hb2 = 1.0,               # 1/s — high-affinity buffer off-rate
        B_L = 500.0,             # µM — total low-affinity buffer
        B_H = 300.0,             # µM — total high-affinity buffer

        # Glutamate release
        alpha_Glu = 1.0,         # dimensionless — max glutamate release
        V_Glu_half = -40.0,      # mV — release half-activation voltage
        V_Glu_slope = 5.0,       # mV — release slope
        tau_Glu = 5.0            # ms — glutamate clearance time constant
    )
end

# ── 2. Initial Conditions ───────────────────────────────────

"""
    rod_dark_state(params)

Return dark-adapted initial conditions for a rod photoreceptor.

# Arguments
- `params`: named tuple from `default_rod_params()`

# Returns
- 19-element state vector corresponding to dark-adapted equilibrium
"""
function rod_dark_state(params)

        # Decompose state vector
        #V, mKv, hKv, mCa, mKCa, Ca_s, Ca_f, Cab_ls, Cab_hs, Cab_lf, Cab_hf,
        #Rh, Rhi, Tr, PDE, Ca_photo, Cab_photo, cGMP, Glu = u

        u0 = zeros(ROD_STATE_VARS)

        # Membrane potential and gating variables
        u0[1] = -36.186 # V
        u0[2] = 0.430 # mKv
        u0[3] = 0.999 # hKv
        u0[4] = 0.436 # mCa
        u0[5] = 0.642 # mKCa

        # Inner segment Ca dynamics
        u0[6] = 0.0966 # Ca_s
        u0[7] = 0.0966 # Ca_f  
        u0[8] = 80.929 # Cab_ls
        u0[9] = 29.068 # Cab_hs
        u0[10] = 80.929 # Cab_lf
        u0[11] = 29.068 # Cab_hf

        # Phototransduction cascade (dark = no activation)
        u0[12] = 0.0 # Rh
        u0[13] = 0.0 # Rhi
        u0[14] = 0.0 # Tr
        u0[15] = 0.0 # PDE

        # Outer segment Ca and cGMP
        u0[16] = 0.3 # Ca_photo
        u0[17] = 34.88 # Cab_photo
        u0[18] = 2.0 # cGMP

        # Glutamate release (compute from steady-state V)
        u0[19] = release_sigmoid(u0[1], params.V_Glu_half, params.V_Glu_slope)

        return u0
end

"""
    cone_dark_state(params)

Return dark-adapted initial conditions for a cone photoreceptor.

# Arguments
- `params`: named tuple from `default_cone_params()`

# Returns
- 6-element state vector corresponding to dark-adapted equilibrium
"""
function cone_dark_state(params)
        # Decompose state vector
        #V, mKv, hKv, mCa, mKCa, Ca_s, Ca_f, Cab_ls, Cab_hs, Cab_lf, Cab_hf,
        #Rh, Rhi, Tr, PDE, Ca_photo, Cab_photo, cGMP, Glu = u

        u0 = zeros(ROD_STATE_VARS)

        # Membrane potential and gating variables
        u0[1] = -36.186 # V
        u0[2] = 0.430 # mKv
        u0[3] = 0.999 # hKv
        u0[4] = 0.436 # mCa
        u0[5] = 0.642 # mKCa

        # Inner segment Ca dynamics
        u0[6] = 0.0966 # Ca_s
        u0[7] = 0.0966 # Ca_f  
        u0[8] = 80.929 # Cab_ls
        u0[9] = 29.068 # Cab_hs
        u0[10] = 80.929 # Cab_lf
        u0[11] = 29.068 # Cab_hf

        # Phototransduction cascade (dark = no activation)
        u0[12] = 0.0 # Rh
        u0[13] = 0.0 # Rhi
        u0[14] = 0.0 # Tr
        u0[15] = 0.0 # PDE

        # Outer segment Ca and cGMP
        u0[16] = 0.3 # Ca_photo
        u0[17] = 34.88 # Cab_photo
        u0[18] = 2.0 # cGMP

        # Glutamate release (compute from steady-state V)
        u0[19] = release_sigmoid(u0[1], params.V_Glu_half, params.V_Glu_slope)

        return u0
end

# ── 3. Auxiliary Functions ──────────────────────────────────

"""
    rate_fraction(scale, x, width, limit_value)

Compute `scale * x / (exp(x/width) - 1)` with numerical safeguard.
Returns `limit_value` when denominator is near zero.
"""
@inline function rate_fraction(scale::Real, x::Real, width::Real, limit_value::Real)
    den = exp(x / width) - 1.0
    if abs(den) < oftype(den, 1.0e-9)
        return oftype(x, limit_value)
    end
    return scale * x / den
end

"""
    release_sigmoid(V, V_half, V_slope)

Steady-state sigmoidal release function for neurotransmitter release.
"""
@inline function release_sigmoid(V::Real, V_half::Real, V_slope::Real)
    return 1.0 / (1.0 + exp(-(V - V_half) / V_slope))
end

# ------------------------------------------------------------------
# 1.  Stimulus and phototransduction helpers
# ------------------------------------------------------------------
Stim(t, t_on, t_off, Φ; hold = 0) = (t_on <= t <= t_off ? Φ : hold)

J∞(g, kg) = g^3 / (g^3 + kg^3)                      # CNG gating
C∞(C, Cae, K) = C > Cae ? (C - Cae)/(C - Cae + K) : 0.0   # clamp at 0

# ------------------------------------------------------------------
# 2.  Voltage‑gated K+ (IKv)
# ------------------------------------------------------------------
αmKv(v) =  5*(100 - v) / ( exp((100 - v)/42) - 1 )
βmKv(v) =  9 * exp(-(v - 20)/40)

αhKv(v) = 0.15 * exp(-v/22)
βhKv(v) = 0.4125 / ( exp((10 - v)/7) + 1 )

# ------------------------------------------------------------------
# 3.  L‑type Ca2+ current (ICa)
# ------------------------------------------------------------------
αmCa(v) =  3*(80 - v) / ( exp((80 - v)/25) - 1 )
βmCa(v) = 10 / ( 1 + exp((v + 38)/7) )

hCa(v)  = exp((40 - v)/18) / ( 1 + exp((40 - v)/18) )

# ------------------------------------------------------------------
# 4.  Ca2+‑activated K+ current (IK(Ca))
# ------------------------------------------------------------------
αmKCa(v) = 15*(80 - v) / ( exp((80 - v)/40) - 1 )
βmKCa(v) = 20 * exp(-v/35)

mKCas(C) = C / (C + 0.3)                             # BK Ca‑gate

# ------------------------------------------------------------------
# 5.  Ca2+‑activated Cl− current (ICl(Ca))
# ------------------------------------------------------------------
mCl(C)  = 1 / ( 1 + exp((0.37 - C)/0.09) )           # note 0.09 !

# ------------------------------------------------------------------
# 6.  Hyperpolarisation‑activated current (Ih) – 5‑state matrix
# ------------------------------------------------------------------
αh(v) =  8  / ( exp((v + 78)/14) + 1 )
βh(v) = 18 / ( exp(-(v + 8)/19) + 1 )                # minus‑sign & 8

hT(v)  = [
   -4αh(v)       βh(v)        0        0        0
    4αh(v)  -(3αh(v)+βh(v))  2βh(v)    0        0
      0        3αh(v)  -(2αh(v)+2βh(v)) 3βh(v)   0
      0          0        2αh(v)  -(αh(v)+3βh(v)) 4βh(v)
      0          0          0        αh(v)    -4βh(v)
]
# ── 4. Mathematical Models ──────────────────────────────────

"""
    photoreceptor_model!(du, u, p, t)

Biophysical rod photoreceptor model with detailed phototransduction cascade,
inner/outer segment Ca dynamics, and multiple ionic currents.

# Arguments
- `du`: derivative vector (19 elements)
- `u`: state vector (19 elements)
- `p`: tuple `(params, Phi, I_feedback)` where:
  - `params`: named tuple from `default_rod_params()`
  - `Phi`: light stimulus (photons/µm²/ms)
  - `I_feedback`: feedback current (pA)
- `t`: time (ms)

# State vector
`u = [V, mKv, hKv, mCa, mKCa, Ca_s, Ca_f, Cab_ls, Cab_hs, Cab_lf, Cab_hf,
      Rh, Rhi, Tr, PDE, Ca_photo, Cab_photo, cGMP, Glu]`
"""
function photoreceptor_model!(du, u, p, t)
    println("t: $t")
    params, Phi, I_feedback = p

    # Decompose state vector
    V, mKv, hKv, mCa, mKCa, Ca_s, Ca_f, Cab_ls, Cab_hs, Cab_lf, Cab_hf,
        Rh, Rhi, Tr, PDE, Ca_photo, Cab_photo, cGMP, Glu = u
    println("Ca_s: $Ca_s")
    # ── Phototransduction cascade ──

    # Light drive: Phi (photons/µm²/ms) → Rh*/s
    Jhv = params.eta * Phi * 1000.0

    T_free = params.T_tot - Tr
    PDE_free = params.PDE_tot - PDE
    cab_photo_free = params.eT - Cab_photo

    # Rhodopsin activation and inactivation (1/s)
    dRh = Jhv - params.alpha1 * Rh + params.alpha2 * Rhi
    dRhi = params.alpha1 * Rh - (params.alpha2 + params.alpha3) * Rhi

    # Transducin and PDE dynamics (1/s)
    dTr = params.epsilon * Rh * T_free - params.beta1 * Tr +
            params.tau2 * PDE - params.tau1 * Tr * PDE_free
    dPDE = params.tau1 * Tr * PDE_free - params.tau2 * PDE

    # Photocurrent
    J = params.J_max * cGMP^3 / (cGMP^3 + 1000.0)
    Iphoto = -J * (1.0 - exp((V - 8.5) / 17.0))

    # Outer segment Ca dynamics (1/s)
    dCa_photo = params.b * J - params.gamma_Ca * (Ca_photo - params.C0) -
                  params.k1 * cab_photo_free * Ca_photo + params.k2 * Cab_photo
    dCab_photo = params.k1 * cab_photo_free * Ca_photo - params.k2 * Cab_photo

    # cGMP synthesis and hydrolysis (1/s)
    cyclase = params.A_max / (1.0 + (Ca_photo / params.K_c)^4)
    dcGMP = cyclase - cGMP * (params.nu + params.sigma * PDE)

    # ── Ionic currents (inner segment) ──
    dmKv = αmKv(V) * (1 - mKv) - βmKv(V) * mKv
    dhKv = αhKv(V) * (1 - hKv) - βhKv(V) * hKv
    # Voltage-gated K+ (delayed rectifier)
    IKv = params.g_Kv * mKv^3 * hKv * (V - params.E_K)
    
    # Voltage-gated Ca2+
    dmCa = αmCa(V) * (1 - mCa) - βmCa(V) * mCa
    dmKCa = αmKCa(V) * (1 - mKCa) - βmKCa(V) * mKCa
    ECa = 12.5 * log(Ca_s / params.Ca_o)
    ICa = params.g_Ca * mCa^4 * hCa(V) * (V - ECa)

    # Ca-activated Cl-
    mCl = 1.0 / (1.0 + exp((0.37 - Ca_s) / 0.09))
    IClCa = params.g_Cl * mCl * (V - params.E_Cl)

    # Ca-activated K+
    mKCa_inf = Ca_s / (Ca_s + 0.3)
    IKCa = params.g_KCa * mKCa^2 * mKCa_inf * (V - params.E_K)

    # Hyperpolarization-activated current
    h_inf = 1.0 / (1.0 + exp((V - params.V_h_half) / params.k_h))
    Ih = params.g_H * h_inf * (V - params.E_H)

    # Leak current
    IL = params.g_L * (V - params.E_L)

    # Ca exchangers
    Iex = params.J_ex_max * Ca_s / (Ca_s + params.K_ex)
    Iex2 = params.J_ex2_max * Ca_s / (Ca_s + params.K_ex2)

    # Membrane voltage derivative (mV/ms)
    dV = (-Iphoto - Ih - IKv - ICa - IClCa - IKCa - IL - Iex - Iex2 + I_feedback) / params.C_m

    # ── Inner segment Ca compartments ──

    cab_ls_free = params.B_L - Cab_ls
    cab_hs_free = params.B_H - Cab_hs
    cab_lf_free = params.B_L - Cab_lf
    cab_hf_free = params.B_H - Cab_hf

    ca_flux = -(ICa + Iex + Iex2) / (2.0 * params.F * params.V1) * 1.0e-6
    diffusion = params.D_Ca * (params.S1 / (params.delta * params.V1)) * (Ca_s - Ca_f)
    # Submembrane Ca (1/s)
    dCa_s = ca_flux - diffusion 
            -params.Lb1 * Ca_s * cab_ls_free + params.Lb2 * Cab_ls
            -params.Hb1 * Ca_s * cab_hs_free + params.Hb2 * Cab_hs
    dCa_f = diffusion 
            - params.Lb1 * Ca_f * cab_lf_free + params.Lb2 * Cab_lf 
            - params.Hb1 * Ca_f * cab_hf_free + params.Hb2 * Cab_hf

    # Ca buffers (1/s)
    dCab_ls = params.Lb1 * Ca_s * cab_ls_free - params.Lb2 * Cab_ls
    dCab_hs = params.Hb1 * Ca_s * cab_hs_free - params.Hb2 * Cab_hs
    dCab_lf = params.Lb1 * Ca_f * cab_lf_free - params.Lb2 * Cab_lf
    dCab_hf = params.Hb1 * Ca_f * cab_hf_free - params.Hb2 * Cab_hf

    # Glutamate release dynamics
    R_glu = release_sigmoid(V, params.V_Glu_half, params.V_Glu_slope)
    dGlu = (params.alpha_Glu * R_glu - Glu) / params.tau_Glu

    du .= [dV, dmKv, dhKv, dmCa, dmKCa, dCa_s, dCa_f, dCab_ls, dCab_hs, dCab_lf, dCab_hf, dRh, dRhi, dTr, dPDE, dCa_photo, dCab_photo, dcGMP, dGlu]

    return nothing
end

# For rod_model!, we need to pass (params, Phi, I_feedback) where Phi varies with time
# So we create a wrapper that computes Phi from the stimulus
function photoreceptor_model_with_stim!(du, u, p, t, stim)
    params, stim, I_feedback = p
    Phi = compute_stimulus(stim, t)
    photoreceptor_model!(du, u, (params, Phi, I_feedback), t)
    return nothing
end

"""
    photoreceptor_K_current(u, params)

Compute total K+ current from a photoreceptor for Müller/RPE K+ sensing.
Supports both rod and cone models via parameter detection.

# Rod version
Computes IKv + IKCa from rod state and parameters.

# Cone version
Computes IKv from cone state and parameters.
"""
function photoreceptor_K_current(u, params::NamedTuple)
    # Detect model type by checking for rod-specific parameters
    if haskey(params, :g_KCa)
        # Rod model - extract relevant state variables
        V = u[ROD_V_INDEX]
        mKv = u[ROD_MKV_INDEX]
        hKv = u[ROD_HKV_INDEX]
        mKCa = u[ROD_MKCA_INDEX]
        Ca_s = u[ROD_CA_S_INDEX]

        IKv = params.g_Kv * mKv^3 * hKv * (V - params.E_K)
        mKCa_inf = Ca_s / (Ca_s + 0.3)
        IKCa = params.g_KCa * mKCa^2 * mKCa_inf * (V - params.E_K)
        return IKv + IKCa
    else
        # Cone model
        V = u[CONE_V_INDEX]
        w_Kv = 1.0 / (1.0 + exp(-(V + 30.0) / 10.0))
        return params.g_Kv * w_Kv * (V - params.E_K)
    end
end
