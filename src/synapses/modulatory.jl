# ============================================================
# modulatory.jl — Dopaminergic modulation
# Spec §4.3
# ============================================================

"""
    dopamine_gain_factor(da_concentration, kappa)

Compute multiplicative gain factor from dopamine.
`kappa` > 0 for facilitation, < 0 for suppression.
Returns: 1 + kappa * [DA]
"""
function dopamine_gain_factor(da::Real, kappa::Real)
    return 1.0 + kappa * da
end

"""
    apply_modulation(g_target, da, kappa)

Return modulated conductance: g' = g * (1 + kappa * [DA]).
"""
function apply_modulation(g_target::Real, da::Real, kappa::Real)
    return g_target * dopamine_gain_factor(da, kappa)
end
