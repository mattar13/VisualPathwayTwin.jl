# ============================================================
# mglur6.jl — mGluR6 sign-inverting synapse
# Spec §3.3.1 — ON-bipolar specific
# ============================================================

"""
    mglur6_conductance(S_mGluR6, params::mGluR6Params)

TRPM1 conductance: high when mGluR6 cascade is inactive (low glutamate / light).
"""
function mglur6_conductance(S_mGluR6::Real, params::mGluR6Params)
    S_clamped = min(max(S_mGluR6, 0.0), 1.0)
    return params.g_TRPM1_max * (1.0 - S_clamped)
end

"""
    mglur6_current(V, S_mGluR6, params::mGluR6Params)

TRPM1 current. Inward (depolarizing) when TRPM1 is open and V < E_TRPM1.
"""
function mglur6_current(V::Real, S_mGluR6::Real, params::mGluR6Params)
    g = mglur6_conductance(S_mGluR6, params)
    return g * (V - params.E_TRPM1)
end

"""
    mglur6_derivative(S, glu_pre, params::mGluR6Params)

Time derivative of the mGluR6 cascade state variable.
"""
function mglur6_derivative(S::Real, glu_pre::Real, params::mGluR6Params)
    return (params.alpha_mGluR6 * glu_pre - S) / params.tau_mGluR6
end
