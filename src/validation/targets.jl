# ============================================================
# targets.jl — Validation target extraction helpers
# ============================================================

"""
    compute_validation_metrics(sim; t_on=sim.col.stimulus.t_on)

Compute coarse validation metrics from a simulation output `sim` (the NamedTuple
returned by `simulate_flash`). Returns a NamedTuple with peak amplitudes and
latencies for a-/b-/c-waves, OP strength, and rod hyperpolarization.
"""
function compute_validation_metrics(sim; t_on::Real=sim.col.stimulus.t_on)
    t = sim.t
    erg = sim.erg
    rod_V = get(sim.cell_voltages, :rod, zeros(length(t), 0))

    function window_indices(t_start::Real, t_end::Real)
        findall(ti -> (ti >= t_start && ti <= t_end), t)
    end

    function extrema_in_window(signal::AbstractVector, t_start::Real, t_end::Real; mode::Symbol)
        idx = window_indices(t_start, t_end)
        if isempty(idx)
            return (value=NaN, time=NaN)
        end
        if mode == :min
            values = signal[idx]
            min_val, rel_idx = findmin(values)
            return (value=min_val, time=t[idx[rel_idx]])
        elseif mode == :max
            values = signal[idx]
            max_val, rel_idx = findmax(values)
            return (value=max_val, time=t[idx[rel_idx]])
        else
            return (value=NaN, time=NaN)
        end
    end

    a_window = (t_on, t_on + 25.0)
    b_window = (t_on + 25.0, t_on + 80.0)
    op_window = (t_on + 25.0, t_on + 150.0)
    c_window = (t_on + 2000.0, min(t_on + 5000.0, t[end]))

    a_peak = extrema_in_window(erg, a_window...; mode=:min)
    b_peak = extrema_in_window(erg, b_window...; mode=:max)
    c_peak = extrema_in_window(erg, c_window...; mode=:max)

    ops = extract_ops(erg, t)
    op_idx = window_indices(op_window...)
    op_amp = isempty(op_idx) ? NaN : maximum(abs.(ops[op_idx]))

    rod_delta = NaN
    if !isempty(rod_V)
        baseline_idx = window_indices(max(t_on - 50.0, t[1]), t_on)
        flash_idx = window_indices(a_window...)
        if !isempty(baseline_idx) && !isempty(flash_idx)
            baseline = mean(rod_V[baseline_idx, :])
            flash_min = minimum(rod_V[flash_idx, :])
            rod_delta = flash_min - baseline
        end
    end

    return (
        a_wave_min=a_peak.value,
        a_wave_latency=a_peak.time - t_on,
        b_wave_max=b_peak.value,
        b_wave_latency=b_peak.time - t_on,
        ops_amplitude=op_amp,
        c_wave_max=c_peak.value,
        c_wave_latency=c_peak.time - t_on,
        rod_hyperpolarization=rod_delta,
    )
end

