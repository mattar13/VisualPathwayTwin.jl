# ============================================================
# parameters.jl — Default parameter sets for all cell types
# Values from spec Tables 10.1, 3.1.4, 4.1
# ============================================================

# --- Photoreceptors (spec §3.1.4) ---

function default_rod_params()
    RodPhotoreceptorParams(
        C_m = 0.02, g_L = 0.35, E_L = -77.0,
        eta = 0.67,
        alpha1 = 50.0, alpha2 = 0.0003, alpha3 = 0.03,
        epsilon = 0.5, beta1 = 2.5, tau1 = 0.2, tau2 = 5.0,
        T_tot = 1000.0, PDE_tot = 100.0,
        J_max = 5040.0, b = 0.25, gamma_Ca = 50.0, C0 = 0.1,
        k1 = 0.2, k2 = 0.8, eT = 500.0,
        A_max = 65.6, K_c = 0.1, nu = 0.4, sigma = 1.0,
        g_Kv = 2.0, E_K = -74.0, g_Ca = 0.7, Ca_o = 1600.0,
        g_Cl = 2.0, E_Cl = -20.0, g_KCa = 5.0,
        alpha_Glu = 1.0, V_Glu_half = -40.0, V_Glu_slope = 5.0, tau_Glu = 5.0
    )
end

function default_cone_params()
    # NOTE: G_dark must equal alpha_G / beta_G for steady state in darkness
    PhototransductionParams(
        C_m = 20.0, g_L = 2.0, E_L = -70.0,
        eta = 0.5, tau_R = 10.0,           # faster inactivation
        alpha_G = 2.5, beta_G = 0.5, gamma_PDE = 5.0,  # FIXED: alpha_G reduced to match G_dark
        n_Ca = 2, Ca_dark = 0.3, G_dark = 5.0,  # weaker Ca feedback
        n_G = 2.5, g_CNG_max = 3.0, E_CNG = 0.0,
        f_Ca = 0.12, k_ex = 0.1, B_Ca = 50.0,
        g_H = 1.0, tau_H = 30.0, V_h_half = -70.0, k_h = -10.0, E_H = -30.0,  # less I_H
        g_Kv = 3.0, E_K = -84.0,
        alpha_Glu = 1.0, V_Glu_half = -40.0, V_Glu_slope = 5.0, tau_Glu = 5.0
    )
end

# --- Morris-Lecar cells (spec Table 10.1) ---

function default_hc_params()
    MLParams(
        C_m = 20.0, g_L = 2.0, g_Ca = 4.0, g_K = 8.0,
        E_L = -60.0, E_Ca = 120.0, E_K = -84.0,
        V1 = -1.2, V2 = 18.0, V3 = 12.0, V4 = 17.0,
        phi = 0.067,
        release = NTReleaseParams(alpha=1.0, V_half=-40.0, V_slope=5.0, tau=10.0)
    )
end

function default_on_bc_params()
    MLParams(
        C_m = 20.0, g_L = 2.0, g_Ca = 4.0, g_K = 8.0,
        E_L = -60.0, E_Ca = 120.0, E_K = -84.0,
        V1 = -1.2, V2 = 18.0, V3 = 12.0, V4 = 17.0,
        phi = 0.067,
        release = NTReleaseParams(alpha=1.0, V_half=-35.0, V_slope=5.0, tau=5.0)
    )
end

function default_off_bc_params()
    MLParams(
        C_m = 20.0, g_L = 2.0, g_Ca = 4.0, g_K = 8.0,
        E_L = -50.0, E_Ca = 120.0, E_K = -84.0,  # E_L = -50 for OFF
        V1 = -1.2, V2 = 18.0, V3 = 2.0, V4 = 17.0,
        phi = 0.067,
        release = NTReleaseParams(alpha=1.0, V_half=-35.0, V_slope=5.0, tau=5.0)
    )
end

function default_a2_params()
    MLParams(
        C_m = 10.0, g_L = 2.0, g_Ca = 8.0, g_K = 12.0,  # fast dynamics
        E_L = -60.0, E_Ca = 120.0, E_K = -84.0,
        V1 = -1.2, V2 = 18.0, V3 = -10.0, V4 = 12.0,  # shifted V3 for oscillations
        phi = 0.2,  # fast
        release = NTReleaseParams(alpha=1.0, V_half=-20.0, V_slope=5.0, tau=4.0)  # glycine
    )
end

function default_gaba_params()
    MLParams(
        C_m = 10.0, g_L = 2.0, g_Ca = 8.0, g_K = 12.0,  # fast dynamics
        E_L = -60.0, E_Ca = 120.0, E_K = -84.0,
        V1 = -1.2, V2 = 18.0, V3 = -8.0, V4 = 12.0,
        phi = 0.15,
        release = NTReleaseParams(alpha=1.0, V_half=-20.0, V_slope=5.0, tau=7.0)  # GABA
    )
end

function default_da_params()
    MLParams(
        C_m = 20.0, g_L = 2.0, g_Ca = 4.0, g_K = 8.0,
        E_L = -60.0, E_Ca = 120.0, E_K = -84.0,
        V1 = -1.2, V2 = 18.0, V3 = 12.0, V4 = 17.0,
        phi = 0.067,
        release = NTReleaseParams(alpha=1.0, V_half=-20.0, V_slope=5.0, tau=200.0)  # slow DA
    )
end

function default_gc_params()
    MLParams(
        C_m = 25.0, g_L = 2.0, g_Ca = 5.0, g_K = 10.0,
        E_L = -65.0, E_Ca = 120.0, E_K = -84.0,
        V1 = -1.2, V2 = 18.0, V3 = 2.0, V4 = 17.0,
        phi = 0.04,
        release = NTReleaseParams()  # GC doesn't release NT in this model
    )
end

# --- Non-neural cells ---

default_muller_params() = MullerParams()
default_rpe_params() = RPEParams()
default_mglur6_params() = mGluR6Params()
default_erg_weights() = ERGWeights()

# --- Column builder ---

function build_retinal_column(; regime::Symbol=:scotopic, kwargs...)
    col = RetinalColumn(
        pop = PopulationSizes(; get(kwargs, :pop, NamedTuple())...),
        rod_params = default_rod_params(),
        cone_params = default_cone_params(),
        hc_params = default_hc_params(),
        on_params = default_on_bc_params(),
        off_params = default_off_bc_params(),
        a2_params = default_a2_params(),
        gaba_params = default_gaba_params(),
        da_params = default_da_params(),
        gc_params = default_gc_params(),
        muller_params = default_muller_params(),
        rpe_params = default_rpe_params(),
        mglur6_params = default_mglur6_params(),
        erg_weights = default_erg_weights(),
        stimulus = StimulusProtocol(; get(kwargs, :stimulus, NamedTuple())...)
    )
    return col
end
