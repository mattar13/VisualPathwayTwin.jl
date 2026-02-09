# ============================================================
# types.jl — Core type definitions for RetinalTwin
# ============================================================

using Parameters

const ROD_STATE_VARS = 19
const CONE_STATE_VARS = 6

const ROD_V_INDEX = 1
const ROD_GLU_INDEX = 19
const CONE_V_INDEX = 4
const CONE_GLU_INDEX = 6

"""
    NTReleaseParams

Parameters for neurotransmitter release dynamics.
Used by cells that release NT (bipolars, amacrines).
"""
@with_kw struct NTReleaseParams
    alpha::Float64 = 1.0        # max release rate
    V_half::Float64 = -35.0     # mV — half-activation voltage
    V_slope::Float64 = 5.0      # mV — sigmoid slope
    tau::Float64 = 5.0          # ms — release/clearance time constant
end

"""
    MLParams

Morris-Lecar parameters for a single cell type.
"""
@with_kw struct MLParams
    C_m::Float64 = 20.0         # pF — membrane capacitance
    g_L::Float64 = 2.0          # nS — leak conductance
    g_Ca::Float64 = 4.0         # nS — calcium conductance
    g_K::Float64 = 8.0          # nS — potassium conductance
    E_L::Float64 = -60.0        # mV — leak reversal
    E_Ca::Float64 = 120.0       # mV — calcium reversal
    E_K::Float64 = -84.0        # mV — potassium reversal
    V1::Float64 = -1.2          # mV — m_inf half-activation
    V2::Float64 = 18.0          # mV — m_inf slope
    V3::Float64 = 2.0           # mV — w_inf half-activation
    V4::Float64 = 30.0          # mV — w_inf slope
    phi::Float64 = 0.04         # dimensionless — w time constant scaling
    release::NTReleaseParams = NTReleaseParams()
end

"""
    PhototransductionParams

Parameters for the 3-compartment phototransduction cascade
plus photoreceptor membrane dynamics.
"""
@with_kw struct PhototransductionParams
    # Membrane
    C_m::Float64 = 20.0         # pF
    g_L::Float64 = 2.0          # nS
    E_L::Float64 = -70.0        # mV
    # Phototransduction cascade
    eta::Float64 = 0.67         # quantum efficiency
    tau_R::Float64 = 80.0       # ms — R* inactivation
    alpha_G::Float64 = 20.0     # µM/ms — basal cGMP synthesis
    beta_G::Float64 = 0.5       # 1/ms — basal cGMP hydrolysis
    gamma_PDE::Float64 = 5.0    # gain of R* on PDE
    n_Ca::Int = 4               # Ca feedback cooperativity
    Ca_dark::Float64 = 0.3      # µM — dark calcium
    G_dark::Float64 = 5.0       # µM — dark cGMP
    n_G::Float64 = 3.0          # Hill coeff for CNG
    g_CNG_max::Float64 = 20.0   # nS — max CNG conductance
    E_CNG::Float64 = 0.0        # mV — CNG reversal
    f_Ca::Float64 = 0.12        # fraction of CNG current carried by Ca
    k_ex::Float64 = 0.1         # 1/ms — Ca extrusion rate
    B_Ca::Float64 = 50.0        # Ca buffering capacity
    # I_H channel
    g_H::Float64 = 2.0          # nS — max I_H conductance
    tau_H::Float64 = 50.0       # ms — I_H time constant
    V_h_half::Float64 = -70.0   # mV — I_H half-activation
    k_h::Float64 = -10.0        # mV — I_H slope (negative: activates with hyperpolarization)
    E_H::Float64 = -30.0        # mV — I_H reversal (mixed Na/K)
    # Voltage-gated K
    g_Kv::Float64 = 3.0         # nS
    E_K::Float64 = -84.0        # mV
    # Glutamate release
    alpha_Glu::Float64 = 1.0    # max glutamate
    V_Glu_half::Float64 = -40.0 # mV
    V_Glu_slope::Float64 = 5.0  # mV
    tau_Glu::Float64 = 5.0      # ms — release time constant
end

"""
    RodPhotoreceptorParams

Parameters for the biophysical rod model in `docs/rod_photoreceptor_model.md`.
Rates are in 1/s and converted to 1/ms inside the update function.
"""
@with_kw struct RodPhotoreceptorParams
    # Membrane
    C_m::Float64 = 0.02          # nF
    g_L::Float64 = 0.35          # nS
    E_L::Float64 = -77.0         # mV

    # Light -> Rh* drive
    eta::Float64 = 0.67

    # Phototransduction cascade
    alpha1::Float64 = 50.0       # 1/s
    alpha2::Float64 = 0.0003     # 1/s
    alpha3::Float64 = 0.03       # 1/s
    epsilon::Float64 = 0.5       # 1/(s*uM)
    beta1::Float64 = 2.5         # 1/s
    tau1::Float64 = 0.2          # 1/(s*uM)
    tau2::Float64 = 5.0          # 1/s
    T_tot::Float64 = 1000.0      # uM
    PDE_tot::Float64 = 100.0     # uM
    J_max::Float64 = 5040.0      # pA
    b::Float64 = 0.25            # uM/(s*pA)
    gamma_Ca::Float64 = 50.0     # 1/s
    C0::Float64 = 0.1            # uM
    k1::Float64 = 0.2            # 1/(s*uM)
    k2::Float64 = 0.8            # 1/s
    eT::Float64 = 500.0          # uM
    A_max::Float64 = 65.6        # uM/s
    K_c::Float64 = 0.1           # uM
    nu::Float64 = 0.4            # 1/s
    sigma::Float64 = 1.0         # 1/(s*uM)

    # Ionic currents
    g_H::Float64 = 1.5           # nS
    E_H::Float64 = -32.0         # mV
    V_h_half::Float64 = -70.0    # mV
    k_h::Float64 = -7.0          # mV

    g_Kv::Float64 = 2.0          # nS
    E_K::Float64 = -74.0         # mV

    g_Ca::Float64 = 0.7          # nS
    Ca_o::Float64 = 1600.0       # uM

    g_Cl::Float64 = 2.0          # nS
    E_Cl::Float64 = -20.0        # mV

    g_KCa::Float64 = 5.0         # nS

    # Exchangers in membrane/Ca equations
    J_ex_max::Float64 = 20.0     # pA
    K_ex::Float64 = 0.2          # uM
    J_ex2_max::Float64 = 5.0     # pA
    K_ex2::Float64 = 0.5         # uM

    # Intracellular calcium system
    F::Float64 = 96485.33212     # C/mol
    V1::Float64 = 3.812e-13      # L
    V2::Float64 = 5.236e-13      # L
    S1::Float64 = 3.142e-8       # effective geometry factor
    delta::Float64 = 0.05
    D_Ca::Float64 = 6.0e-8

    Lb1::Float64 = 2.0           # 1/(s*uM)
    Lb2::Float64 = 1.0           # 1/s
    Hb1::Float64 = 1.11          # 1/(s*uM)
    Hb2::Float64 = 1.0           # 1/s
    B_L::Float64 = 500.0         # uM
    B_H::Float64 = 300.0         # uM

    # Glutamate release state for downstream circuitry
    alpha_Glu::Float64 = 1.0
    V_Glu_half::Float64 = -40.0
    V_Glu_slope::Float64 = 5.0
    tau_Glu::Float64 = 5.0       # ms
end

"""
    SynapseParams

Parameters for a single synaptic connection type.
"""
@with_kw struct SynapseParams
    g_max::Float64              # nS — maximal conductance
    E_rev::Float64              # mV — reversal potential (NaN for modulatory)
    tau_s::Float64              # ms — postsynaptic gating time constant
    nt_type::Symbol = :E        # :E, :I, or :M
    receptor::Symbol = :iGluR   # :iGluR, :mGluR6, :GlyR, :GABA_A, :D1R, :feedback
end

"""
    mGluR6Params

Parameters specific to the mGluR6 sign-inverting synapse.
"""
@with_kw struct mGluR6Params
    g_TRPM1_max::Float64 = 10.0  # nS
    E_TRPM1::Float64 = 0.0       # mV (non-selective cation)
    alpha_mGluR6::Float64 = 1.0   # scaling
    tau_mGluR6::Float64 = 30.0    # ms — metabotropic time constant
end

"""
    MullerParams

Müller glial cell parameters.
"""
@with_kw struct MullerParams
    C_m::Float64 = 30.0          # pF
    g_Kir_end::Float64 = 5.0     # nS — endfoot Kir conductance
    g_Kir_stalk::Float64 = 2.0   # nS — stalk Kir conductance
    K_o_rest::Float64 = 3.0      # mM — resting extracellular K+
    K_i::Float64 = 140.0         # mM — intracellular K+
    tau_K_diffusion::Float64 = 200.0 # ms — K+ clearance
    alpha_K::Float64 = 0.001     # K+ release per unit current
    V_max_EAAT::Float64 = 0.5    # µM/ms
    K_m_EAAT::Float64 = 10.0     # µM
end

"""
    RPEParams

RPE parameters for c-wave generation.
"""
@with_kw struct RPEParams
    tau_RPE::Float64 = 3000.0    # ms — very slow dynamics
    g_K_apical::Float64 = 5.0    # nS
    g_Cl_baso::Float64 = 2.0     # nS
    g_L_RPE::Float64 = 0.5       # nS
    E_Cl::Float64 = -50.0        # mV
    E_L_RPE::Float64 = -60.0     # mV
    K_sub_rest::Float64 = 3.0    # mM
    k_RPE::Float64 = 0.0005      # K+ clearance rate
    alpha_K_RPE::Float64 = 0.001 # K+ flux scaling
    K_i::Float64 = 140.0         # mM — intracellular K+
end

"""
    ERGWeights

Weights for computing the ERG field potential from cell currents.
"""
@with_kw struct ERGWeights
    rod::Float64 = 1.0
    cone::Float64 = 0.5
    on_bc::Float64 = -2.0
    off_bc::Float64 = -1.0
    a2::Float64 = 0.3
    gaba::Float64 = 0.3
    da::Float64 = 0.05
    gc::Float64 = 0.1
    muller::Float64 = -1.5
    rpe::Float64 = -1.0
end

"""
    StimulusProtocol

Light stimulus specification.
"""
@with_kw struct StimulusProtocol
    I_0::Float64 = 1000.0       # photons/µm²/ms
    t_on::Float64 = 200.0       # ms — flash onset
    t_dur::Float64 = 10.0       # ms — flash duration
    background::Float64 = 0.0   # background illumination
end

"""
    PopulationSizes

Cell population sizes for the retinal column.
"""
@with_kw struct PopulationSizes
    n_rod::Int = 20
    n_cone::Int = 5
    n_hc::Int = 2
    n_on::Int = 1
    n_off::Int = 1
    n_a2::Int = 3
    n_gaba::Int = 3
    n_dopa::Int = 1
    n_gc::Int = 1
    n_muller::Int = 1
    n_rpe::Int = 1
end

"""
    RetinalColumn

Complete parameters for one retinal column.
"""
@with_kw struct RetinalColumn
    pop::PopulationSizes = PopulationSizes()
    rod_params::RodPhotoreceptorParams = RodPhotoreceptorParams()
    cone_params::PhototransductionParams = PhototransductionParams()
    hc_params::MLParams = MLParams()
    on_params::MLParams = MLParams()
    off_params::MLParams = MLParams()
    a2_params::MLParams = MLParams()
    gaba_params::MLParams = MLParams()
    da_params::MLParams = MLParams()
    gc_params::MLParams = MLParams()
    muller_params::MullerParams = MullerParams()
    rpe_params::RPEParams = RPEParams()
    mglur6_params::mGluR6Params = mGluR6Params()
    erg_weights::ERGWeights = ERGWeights()
    stimulus::StimulusProtocol = StimulusProtocol()
end

"""
    StateIndex

Named indices into the flat state vector.
Computed from population sizes at construction time.
"""
struct StateIndex
    rod::UnitRange{Int}       # [V, mKv, hKv, mCa, mKCa, Ca_s, Ca_f, Cab_ls, Cab_hs, Cab_lf, Cab_hf, Rh, Rhi, Tr, PDE, Ca_photo, Cab_photo, cGMP, Glu] × N_rod
    cone::UnitRange{Int}      # [R*, G, Ca, V, h, Glu] × N_cone
    hc::UnitRange{Int}        # [V, w, s_Glu] × N_hc
    on_bc::UnitRange{Int}     # [V, w, S_mGluR6, Glu] × N_on
    off_bc::UnitRange{Int}    # [V, w, s_Glu, Glu_release] × N_off
    a2::UnitRange{Int}        # [V, w, Gly] × N_a2
    gaba_ac::UnitRange{Int}   # [V, w, GABA] × N_gaba
    da_ac::UnitRange{Int}     # [V, w, DA] × N_dopa
    gc::UnitRange{Int}        # [V, w] × N_gc
    muller::UnitRange{Int}    # [V_M, K_o_end, K_o_stalk, Glu_o] × N_muller
    rpe::UnitRange{Int}       # [V_RPE, K_sub] × N_rpe
    total::Int
end

function StateIndex(col::RetinalColumn)
    p = col.pop
    idx = Ref(1)

    function next_range(n_cells, vars_per_cell)
        start = idx[]
        idx[] += n_cells * vars_per_cell
        return start:(idx[] - 1)
    end

    rod     = next_range(p.n_rod, ROD_STATE_VARS)
    cone    = next_range(p.n_cone, CONE_STATE_VARS)
    hc      = next_range(p.n_hc, 3)
    on_bc   = next_range(p.n_on, 4)
    off_bc  = next_range(p.n_off, 4)
    a2      = next_range(p.n_a2, 3)
    gaba_ac = next_range(p.n_gaba, 3)
    da_ac   = next_range(p.n_dopa, 3)
    gc      = next_range(p.n_gc, 2)
    muller  = next_range(p.n_muller, 4)
    rpe     = next_range(p.n_rpe, 2)

    return StateIndex(rod, cone, hc, on_bc, off_bc, a2, gaba_ac, da_ac, gc, muller, rpe, idx[] - 1)
end
