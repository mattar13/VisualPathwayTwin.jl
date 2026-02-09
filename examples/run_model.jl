# ============================================================
# run_model.jl — Step-by-step retinal column simulation
#
# Everything is explicit here. No hidden helper functions.
# Run section by section in the REPL to see exactly what's happening.
# ============================================================
using Revise
using RetinalTwin
#If we have revised the code, 
using DifferentialEquations
using Setfield
using GLMakie

function plot_erg_gl(result; show_components::Bool=true, time_range=nothing)
    t = result.t
    erg = result.erg
    comps = result.erg_components

    if time_range !== nothing
        mask = (t .>= time_range[1]) .& (t .<= time_range[2])
        t = t[mask]
        erg = erg[mask]
        comps = Dict(k => v[mask] for (k, v) in comps)
    end

    if show_components
        fig = Figure(size=(900, 800))
        ax1 = Axis(fig[1, 1], title="Total ERG", xlabel="Time (ms)", ylabel="Amplitude (a.u.)")
        lines!(ax1, t, erg, color=:black, linewidth=2)
        hlines!(ax1, [0.0], color=:gray, linestyle=:dash)

        ax2 = Axis(fig[2, 1], title="ERG Components", xlabel="Time (ms)", ylabel="Amplitude (a.u.)")
        component_colors = Dict(
            :a_wave => :blue, :b_wave => :red, :d_wave => :orange,
            :OPs => :green, :P3 => :purple, :c_wave => :brown, :ganglion => :gray
        )
        component_labels = Dict(
            :a_wave => "a-wave (PR)", :b_wave => "b-wave (ON-BC)", :d_wave => "d-wave (OFF-BC)",
            :OPs => "OPs (Amacrine)", :P3 => "P3 (Müller)", :c_wave => "c-wave (RPE)",
            :ganglion => "GC"
        )

        for (key, trace) in comps
            col = get(component_colors, key, :black)
            lab = get(component_labels, key, String(key))
            lines!(ax2, t, trace, color=col, label=lab)
        end
        axislegend(ax2, position=:rt)
        hlines!(ax2, [0.0], color=:gray, linestyle=:dash)
    else
        fig = Figure(size=(900, 400))
        ax = Axis(fig[1, 1], title="ERG", xlabel="Time (ms)", ylabel="Amplitude (a.u.)")
        lines!(ax, t, erg, color=:black, linewidth=2)
        hlines!(ax, [0.0], color=:gray, linestyle=:dash)
    end

    return fig
end

function plot_cell_voltages_gl(result; cell_types::Vector{Symbol}=[:rod, :on_bc, :a2, :gc],
                               time_range=nothing)
    t = result.t
    voltages = result.cell_voltages

    if time_range !== nothing
        mask = (t .>= time_range[1]) .& (t .<= time_range[2])
        t = t[mask]
    else
        mask = trues(length(t))
    end

    n_panels = length(cell_types)
    fig = Figure(size=(900, 250 * n_panels))

    cell_labels = Dict(
        :rod => "Rod", :cone => "Cone", :hc => "Horizontal Cell",
        :on_bc => "ON-Bipolar", :off_bc => "OFF-Bipolar",
        :a2 => "A2 Amacrine", :gaba_ac => "GABA Amacrine",
        :da_ac => "DA Amacrine", :gc => "Ganglion Cell",
        :muller => "Müller Glia", :rpe => "RPE"
    )

    for (idx, ct) in enumerate(cell_types)
        label = get(cell_labels, ct, String(ct))
        ax = Axis(fig[idx, 1], title=label, xlabel="Time (ms)", ylabel="V (mV)")

        if haskey(voltages, ct)
            V = voltages[ct][mask, :]
            n_cells = size(V, 2)
            for ci in 1:min(n_cells, 5)
                lines!(ax, t, V[:, ci], linewidth=1)
            end
            if n_cells > 5
                lines!(ax, t, vec(mean(V, dims=2)), color=:black, linewidth=2, linestyle=:dash)
            end
        end
    end

    return fig
end

function plot_phototransduction_breakdown(result; cell::Symbol=:rod, cell_index::Int=1, time_range=nothing)
    t = result.t
    sol = result.solution
    stimulus = [compute_stimulus(result.col.stimulus, ti) for ti in t]

    if cell == :rod
        base = result.sidx.rod[1] + (cell_index - 1) * 6
    elseif cell == :cone
        base = result.sidx.cone[1] + (cell_index - 1) * 6
    else
        error("Phototransduction breakdown only supports :rod or :cone, got $(cell).")
    end

    r_star = [u[base] for u in sol.u]
    g_state = [u[base + 1] for u in sol.u]
    ca = [u[base + 2] for u in sol.u]
    v = [u[base + 3] for u in sol.u]
    h = [u[base + 4] for u in sol.u]
    glu = [u[base + 5] for u in sol.u]

    if time_range !== nothing
        mask = (t .>= time_range[1]) .& (t .<= time_range[2])
        t = t[mask]
        stimulus = stimulus[mask]
        r_star = r_star[mask]
        g_state = g_state[mask]
        ca = ca[mask]
        v = v[mask]
        h = h[mask]
        glu = glu[mask]
    end

    fig = Figure(size=(950, 900))
    axes = Axis[]

    push!(axes, Axis(fig[1, 1], title="Phototransduction breakdown ($(uppercase(String(cell))) $(cell_index))",
                     xlabel="Time (ms)", ylabel="Stimulus (ph/μm²/ms)"))
    lines!(axes[end], t, stimulus, color=:orange, linewidth=2)

    push!(axes, Axis(fig[2, 1], xlabel="Time (ms)", ylabel="R* (a.u.)"))
    lines!(axes[end], t, r_star, color=:blue, linewidth=2)

    push!(axes, Axis(fig[3, 1], xlabel="Time (ms)", ylabel="G (a.u.)"))
    lines!(axes[end], t, g_state, color=:purple, linewidth=2)

    push!(axes, Axis(fig[4, 1], xlabel="Time (ms)", ylabel="Ca (a.u.)"))
    lines!(axes[end], t, ca, color=:teal, linewidth=2)

    push!(axes, Axis(fig[5, 1], xlabel="Time (ms)", ylabel="V (mV)"))
    lines!(axes[end], t, v, color=:black, linewidth=2)

    push!(axes, Axis(fig[6, 1], xlabel="Time (ms)", ylabel="Glu (a.u.)"))
    lines!(axes[end], t, glu, color=:red, linewidth=2)

    linkxaxes!(axes...)

    return fig
end

# ── 1. Build the retinal column ──────────────────────────────

col = build_retinal_column(regime=:scotopic)
sidx = StateIndex(col)

println("Retinal column: $(sidx.total) ODEs")
println("  Rods: $(col.pop.n_rod)  |  Cones: $(col.pop.n_cone)")
println("  HC: $(col.pop.n_hc)  |  ON-BC: $(col.pop.n_on)  |  OFF-BC: $(col.pop.n_off)")
println("  A2: $(col.pop.n_a2)  |  GABA: $(col.pop.n_gaba)  |  DA: $(col.pop.n_dopa)")
println("  GC: $(col.pop.n_gc)  |  Müller: $(col.pop.n_muller)  |  RPE: $(col.pop.n_rpe)")

# State vector layout (flat vector, 1-indexed ranges):
println("\nState index ranges:")
println("  rod:     $(sidx.rod)      (6 vars × $(col.pop.n_rod): R*, G, Ca, V, h, Glu)")
println("  cone:    $(sidx.cone)  (6 vars × $(col.pop.n_cone): R*, G, Ca, V, h, Glu)")
println("  hc:      $(sidx.hc)  (3 vars × $(col.pop.n_hc): V, w, s_Glu)")
println("  on_bc:   $(sidx.on_bc)  (4 vars × $(col.pop.n_on): V, w, S_mGluR6, Glu)")
println("  off_bc:  $(sidx.off_bc)  (4 vars × $(col.pop.n_off): V, w, s_Glu, Glu)")
println("  a2:      $(sidx.a2)  (3 vars × $(col.pop.n_a2): V, w, Gly)")
println("  gaba_ac: $(sidx.gaba_ac)  (3 vars × $(col.pop.n_gaba): V, w, GABA)")
println("  da_ac:   $(sidx.da_ac)  (3 vars × $(col.pop.n_dopa): V, w, DA)")
println("  gc:      $(sidx.gc)  (2 vars × $(col.pop.n_gc): V, w)")
println("  muller:  $(sidx.muller)  (4 vars × $(col.pop.n_muller): V_M, K_o_end, K_o_stalk, Glu_o)")
println("  rpe:     $(sidx.rpe)  (2 vars × $(col.pop.n_rpe): V_RPE, K_sub)")

# ── 2. Configure stimulus ────────────────────────────────────

intensity = 10.0      # photons/µm²/ms (dim scotopic flash)
duration  = 10.0      # ms
t_on      = 200.0     # ms — flash onset (give 200 ms dark baseline)

col = @set col.stimulus.I_0    = intensity
col = @set col.stimulus.t_dur  = duration
col = @set col.stimulus.t_on   = t_on

println("\nStimulus: $(intensity) ph/µm²/ms, onset=$(t_on) ms, duration=$(duration) ms")
println("  Background: $(col.stimulus.background)")

# Verify the stimulus function works:
println("  Phi at t=0:   $(compute_stimulus(col.stimulus, 0.0))")
println("  Phi at t=205: $(compute_stimulus(col.stimulus, 205.0))")
println("  Phi at t=220: $(compute_stimulus(col.stimulus, 220.0))")

# ── 3. Connection table ──────────────────────────────────────

conns = default_connections()
println("\n$(length(conns)) synaptic connections:")
for c in conns
    println("  $(c.pre) → $(c.post)  [$(c.nt_type), $(c.receptor), g=$(c.g_max) nS, E=$(c.E_rev) mV]")
end

# ── 4. Initial conditions (dark-adapted) ─────────────────────

u0 = dark_adapted_state(col, sidx)
println("\nInitial state vector: $(length(u0)) elements")

# Spot-check a few values:
rod1_offset = sidx.rod[1]
println("  Rod 1: R*=$(u0[rod1_offset]), G=$(u0[rod1_offset+1]), Ca=$(u0[rod1_offset+2]), V=$(u0[rod1_offset+3]) mV, h=$(u0[rod1_offset+4]), Glu=$(u0[rod1_offset+5])")

on_offset = sidx.on_bc[1]
println("  ON-BC: V=$(u0[on_offset]) mV, w=$(u0[on_offset+1]), S_mGluR6=$(u0[on_offset+2]), Glu=$(u0[on_offset+3])")

gc_offset = sidx.gc[1]
println("  GC:    V=$(u0[gc_offset]) mV, w=$(u0[gc_offset+1])")

println("  All finite? $(all(isfinite, u0))")
println("  Any NaN?    $(any(isnan, u0))")

# ── 5. Set up and solve the ODE ──────────────────────────────

t_total  = 1000.0   # ms simulation window
dt_save  = 0.1      # ms save interval
tspan    = (0.0, t_total)
p        = (col, sidx, conns)

println("\nSetting up ODE problem...")
println("  tspan: $(tspan)")
println("  dt_save: $(dt_save) ms → expect ~$(Int(t_total / dt_save)) saved timepoints")

prob = ODEProblem(retinal_column_rhs!, u0, tspan, p)

# Test the RHS at t=0 to make sure it doesn't error or produce NaN:
du_test = similar(u0)
retinal_column_rhs!(du_test, u0, p, 0.0)
println("  RHS at t=0: all finite? $(all(isfinite, du_test))")
println("  RHS at t=0: any NaN?    $(any(isnan, du_test))")
println("  RHS max |du|: $(maximum(abs, du_test))")

# Also test at flash onset:
retinal_column_rhs!(du_test, u0, p, t_on + 1.0)
println("  RHS at t=$(t_on+1): all finite? $(all(isfinite, du_test))")
println("  RHS at t=$(t_on+1): max |du|: $(maximum(abs, du_test))")

println("\nSolving ODE (Rodas5 — stiff implicit solver)...")
sol = solve(prob, Rodas5();
            saveat=dt_save,
            abstol=1e-6, reltol=1e-4,
            maxiters=1_000_000)

println("  Solver return code: $(sol.retcode)")
println("  Timepoints saved:  $(length(sol.t))")
println("  Time range:        $(sol.t[1]) → $(sol.t[end]) ms")
println("  Solution finite?   $(all(u -> all(isfinite, u), sol.u))")

if length(sol.t) < 10
    println("\n⚠ WARNING: Very few timepoints returned. The solver may have failed.")
    println("  sol.t = $(sol.t)")
    if length(sol.u) >= 1
        u_last = sol.u[end]
        println("  Last state has NaN? $(any(isnan, u_last))")
        println("  Last state has Inf? $(any(isinf, u_last))")
        # Find which variables went bad:
        bad_idx = findall(x -> !isfinite(x), u_last)
        if !isempty(bad_idx)
            println("  Non-finite at indices: $(bad_idx)")
        end
    end
    println("\nStopping here — fix the solver issue before continuing.")
else
    # ── 6. Inspect raw solution ──────────────────────────────────

    println("\nSpot-checking solution at key timepoints:")
    for t_check in [0.0, 100.0, t_on, t_on + 5.0, t_on + 50.0, t_on + 200.0, 500.0]
        ti = argmin(abs.(sol.t .- t_check))
        u = sol.u[ti]
        rod_V = u[rod1_offset + 3]
        on_V  = u[on_offset]
        gc_V  = u[gc_offset]
        println("  t=$(round(sol.t[ti], digits=1)) ms:  Rod V=$(round(rod_V, digits=2)) mV,  ON-BC V=$(round(on_V, digits=2)) mV,  GC V=$(round(gc_V, digits=2)) mV")
    end

    # ── 7. Compute ERG ───────────────────────────────────────────

    println("\nComputing ERG field potential...")
    erg, components = compute_erg(sol, col, sidx)
    println("  ERG trace length: $(length(erg))")
    println("  Components: $(keys(components))")
    println("  a-wave min:  $(round(minimum(components[:a_wave]), digits=4))")
    println("  b-wave max:  $(round(maximum(components[:b_wave]), digits=4))")
    println("  OP range:    $(round(minimum(components[:OPs]), digits=4)) to $(round(maximum(components[:OPs]), digits=4))")
    println("  c-wave max:  $(round(maximum(components[:c_wave]), digits=4))")

    # ── 8. Extract voltages ──────────────────────────────────────

    println("\nExtracting cell voltages...")
    voltages = extract_voltages(sol, col, sidx)
    for (name, V) in voltages
        println("  $(name): $(size(V, 2)) cells, V range = $(round(minimum(V), digits=1)) to $(round(maximum(V), digits=1)) mV")
    end

    # ── 9. Build result tuple (for plotting functions) ───────────

    result = (t=sol.t, erg=erg, erg_components=components,
              cell_voltages=voltages, solution=sol, col=col, sidx=sidx)

    # ── 10. Plot with GLMakie ────────────────────────────────────
    show_plots = true
    if show_plots
        fig_erg = plot_erg_gl(result, show_components=true)
        display(fig_erg)

        fig_v = plot_cell_voltages_gl(result,
            cell_types=[:rod, :cone, :on_bc, :off_bc, :a2, :gaba_ac, :gc])
        display(fig_v)

        fig_photo = plot_phototransduction_breakdown(result, cell=:rod, cell_index=1)
        display(fig_photo)
    end

    println("\nDone. Result tuple available as `result`.")
    println("GLMakie plots displayed: ERG, cell voltages, phototransduction breakdown")
end
