# ============================================================
# synapse.jl — Synaptic transmission utilities
# Spec §2.3, §2.4, §4.2
# ============================================================

"""
    nt_release_sigmoid(V, V_half, V_slope)

Sigmoid activation for neurotransmitter release: T_∞(V).
"""
function nt_release_sigmoid(V::Real, V_half::Real, V_slope::Real)
    return 1.0 / (1.0 + exp(-(V - V_half) / V_slope))
end

"""
    mean_nt(u, range, nt_var_offset, vars_per_cell, n_cells)

Compute population-averaged neurotransmitter concentration from flat state vector.
`nt_var_offset` is the 1-based index of the NT variable within a single cell's state.
"""
function mean_nt(u, range::UnitRange{Int}, nt_var_offset::Int,
                 vars_per_cell::Int, n_cells::Int)
    if n_cells == 0
        return 0.0
    end
    total = 0.0
    base = range[1]
    for i in 0:(n_cells - 1)
        total += u[base + i * vars_per_cell + (nt_var_offset - 1)]
    end
    return total / n_cells
end

"""
    weighted_mean(a, n_a, b, n_b)

Weighted average of two population NT concentrations.
"""
function weighted_mean(a::Real, n_a::Int, b::Real, n_b::Int)
    total_n = n_a + n_b
    if total_n == 0
        return 0.0
    end
    return (a * n_a + b * n_b) / total_n
end

"""
    synaptic_current(g_max, s, V, E_rev)

Standard ionotropic synaptic current: I = g_max * s * (V - E_rev).
"""
function synaptic_current(g_max::Real, s::Real, V::Real, E_rev::Real)
    return g_max * s * (V - E_rev)
end

"""
    mean_voltage(u, range, v_offset, vars_per_cell, n_cells)

Compute population-averaged voltage from flat state vector.
"""
function mean_voltage(u, range::UnitRange{Int}, v_offset::Int,
                      vars_per_cell::Int, n_cells::Int)
    if n_cells == 0
        return -60.0  # default resting
    end
    total = 0.0
    base = range[1]
    for i in 0:(n_cells - 1)
        total += u[base + i * vars_per_cell + (v_offset - 1)]
    end
    return total / n_cells
end
