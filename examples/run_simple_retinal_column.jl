# ============================================================
# simple_retinal_column.jl
# Simple example using retinal_column_model! for photoreceptor → ON bipolar coupling
# ============================================================
using Revise
using RetinalTwin
using DifferentialEquations
using CairoMakie

println("=" ^ 60)
println("Simple Retinal Column: Photoreceptor → ON Bipolar")
println("=" ^ 60)

# ── 1. Load parameters ───────────────────────────────────────
println("\n1. Loading parameters...")
params = default_retinal_params()
println("  ✓ Photoreceptor kHYDRO: $(params.PHOTORECEPTOR_PARAMS.kHYDRO)")
println("  ✓ Photoreceptor parameters: $(length(fieldnames(typeof(params.PHOTORECEPTOR_PARAMS)))) fields")
println("  ✓ ON bipolar parameters: $(length(fieldnames(typeof(params.ON_BIPOLAR_PARAMS)))) fields")

params.PHOTORECEPTOR_PARAMS
# ── 2. Build initial conditions ──────────────────────────────
println("\n2. Building initial conditions...")
u0 = retinal_column_initial_conditions(params)

println("  ✓ Total state variables: $(length(u0))")

# ── 3. Define stimulus ────────────────────────────────────────
println("\n3. Computing steady state (dark, no stimulus)...")
#Create an empty stimulus function
stim_func_dark(t) = RetinalTwin.single_flash(t;photon_flux = 0.0)

tspan = (0.0, 1000.0)
prob = ODEProblem(retinal_column_model!, u0, tspan, (params, stim_func_dark))
sol = solve(prob, Rodas5(); save_everystep=false, save_start=false,save_end=true, abstol=1e-6, reltol=1e-4)
u0 = sol.u[end]

stim_start = 50.0   # ms
stim_end = 55.0      # ms (5 ms flash)
photon_flux = 0.5   # photons/µm²/ms
stim_func(t) = RetinalTwin.single_flash(t; stim_start = stim_start, stim_end = stim_end, photon_flux = photon_flux)
println("\n4. Stimulus configuration:")
println("  ✓ Intensity: $(photon_flux) ph/µm²/ms")
println("  ✓ Duration: $(stim_start) - $(stim_end) ms")

# ── 4. Solve ODE system ───────────────────────────────────────
println("\n5. Solving coupled system...")

tspan = (0.0, 2000.0)  # ms
prob = ODEProblem(retinal_column_model!, u0, tspan, (params, stim_func))

sol = solve(prob, Rodas5(); tstops = [stim_start, stim_end], saveat=1.0, abstol=1e-6, reltol=1e-4)

println("  ✓ Solver status: $(sol.retcode)")
println("  ✓ Time points: $(length(sol.t))")

# ── 5. Extract results ────────────────────────────────────────
println("\n6. Extracting results...")

t = sol.t
t_series = range(0.0, sol.t[end], length=1000)

# Photoreceptor variables
V_photo = map(t-> sol(t)[20], t_series)
Ca_photo = map(t-> sol(t)[15], t_series)
Glu_photo = map(t-> sol(t)[21], t_series)  # Glutamate is now a state variable!

# ON bipolar variables
V_onbc = map(t-> sol(t)[22], t_series)  # V
Ca_onbc = map(t-> sol(t)[24], t_series)  # Ca
Glu_onbc = map(t-> sol(t)[25], t_series)  # ON BC glutamate release

V_offbc = map(t-> sol(t)[28], t_series)  # OFF BC V
Ca_offbc = map(t-> sol(t)[31], t_series)  # OFF BC Ca
A_offbc = map(t-> sol(t)[32], t_series)  # OFF BC iGluR activation
D_offbc = map(t-> sol(t)[33], t_series)  # OFF BC desensitization
Glu_offbc = map(t-> sol(t)[34], t_series)  # OFF BC glutamate release

V_a2 = map(t-> sol(t)[35], t_series)  # A2 amacrine V
Ca_a2 = map(t-> sol(t)[38], t_series)  # A2 amacrine Ca
Gly_a2 = map(t-> sol(t)[41], t_series)  # A2 amacrine glycine release

#%% ── 6. Visualization ──────────────────────────────────────────
println("\n7. Creating plots...")

fig1 = Figure(size=(1200, 500))

# Photoreceptor voltage
ax1a = Axis(fig1[1, 1],
           xlabel="Time (ms)",
           ylabel="V (mV)",
           title="Photoreceptor Membrane Potential")
lines!(ax1a, t_series, V_photo, color=:blue, linewidth=2, label="Photoreceptor")
vlines!(ax1a, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)

ax2a = Axis(fig1[2, 1],
           xlabel="Time (ms)",
           ylabel="Ca (µM)",
           title="Photoreceptor Calcium")
lines!(ax2a, t_series, Ca_photo, color=:purple, linewidth=2)
vlines!(ax2a, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)

# Glutamate release
ax3a = Axis(fig1[3, 1],
           xlabel="Time (ms)",
           ylabel="Glutamate (µM)",
           title="Glutamate Release")
lines!(ax3a, t_series, Glu_photo, color=:orange, linewidth=2)
vlines!(ax3a, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)

ax1b = Axis(fig1[1, 2],
           xlabel="Time (ms)",
           ylabel="V (mV)",
           title="ON Bipolar Membrane Potential")
lines!(ax1b, t_series, V_onbc, color=:green, linewidth=2)
vlines!(ax1b, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)


ax2b = Axis(fig1[2, 2],
           xlabel="Time (ms)",
           ylabel="Ca (µM)",
           title="ON Bipolar Calcium")
lines!(ax2b, t_series, Ca_onbc, color=:purple, linewidth=2)
vlines!(ax2b, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)

ax3b = Axis(fig1[3, 2],
           xlabel="Time (ms)",
           ylabel="Glutamate (µM)",
           title="Glutamate Release")
lines!(ax3b, t_series, Glu_onbc, color=:green, linewidth=2)
vlines!(ax3b, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)


ax1c = Axis(fig1[1, 3],
           xlabel="Time (ms)",
           ylabel="V (mV)",
           title="OFF Bipolar Membrane Potential")
lines!(ax1c, t_series, V_offbc, color=:red, linewidth=2)
vlines!(ax1c, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)

ax2c = Axis(fig1[2, 3],
           xlabel="Time (ms)",
           ylabel="Ca (µM)",
           title="OFF Bipolar Calcium")
lines!(ax2c, t_series, Ca_offbc, color=:purple, linewidth=2)
vlines!(ax2c, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)

ax3c = Axis(fig1[3, 3],
           xlabel="Time (ms)",
           ylabel="Glutamate (µM)",
           title="Glutamate Release")
lines!(ax3c, t_series, Glu_offbc, color=:red, linewidth=2)
vlines!(ax3c, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)

ax1d = Axis(fig1[1, 4],
           xlabel="Time (ms)",
           ylabel="V (mV)",
           title="A2 Amacrine Membrane Potential")
lines!(ax1d, t_series, V_a2, color=:blue, linewidth=2)
vlines!(ax1d, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)

ax2d = Axis(fig1[2, 4],
           xlabel="Time (ms)",
           ylabel="Ca (µM)",
           title="A2 Amacrine Calcium")
lines!(ax2d, t_series, Ca_a2, color=:purple, linewidth=2)
vlines!(ax2d, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)

ax3d = Axis(fig1[3, 4],
           xlabel="Time (ms)",
           ylabel="Glutamate (µM)",
           title="A2 Amacrine Gly release")
lines!(ax3d, t_series, Gly_a2, color=:orange, linewidth=2)
        vlines!(ax3d, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)

save("examples/plots/simple_retinal_column.png", fig1)


#Voltage vs glutamate release
C_rng = range(0.0, 1.0, length=100)
G_rng = [RetinalTwin.R_inf(c, params.ON_BIPOLAR_PARAMS.K_Release, params.ON_BIPOLAR_PARAMS.n_Release) for c in C_rng]
fig2 = Figure(size=(500, 500))
ax3a = Axis(fig2[1, 1],
           xlabel="Ca (µM)",
           ylabel="Glutamate Release (µM)",
           title="Glutamate Release vs Ca")
lines!(ax3a, C_rng, G_rng, color=:blue, linewidth=2)
save("examples/plots/voltage_vs_glutamate_release.png", fig2)


#Plot the desensitization of the OFF bipolar iGluR
fig3 = Figure(size=(500, 500))
ax1a = Axis(fig3[1, 1],
           xlabel="Time (ms)",
           ylabel="V (mV)",
           title="OFF Bipolar Membrane Potential")
lines!(ax1a, t_series, V_offbc, color=:black, linewidth=2)
vlines!(ax1a, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)

ax2a = Axis(fig3[2, 1],
           xlabel="Time (ms)",
           ylabel="Activation",
           title="AMPA/KAR-like activation")
lines!(ax2a, t_series, A_offbc, color=:green, linewidth=2)
vlines!(ax2a, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)


ax3a = Axis(fig3[3, 1],
           xlabel="Time (ms)",
           ylabel="Desensitization",
           title="AMPA/KAR-like desensitization")
lines!(ax3a, t_series, D_offbc, color=:red, linewidth=2)
vlines!(ax3a, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)


ax4a = Axis(fig3[4, 1],
           xlabel="Time (ms)",
           ylabel="Glutamate Release (µM)",
           title="OFF Bipolar Glutamate Release")
lines!(ax4a, t_series, A_offbc .* D_offbc, color=:orange, linewidth=2)
vlines!(ax4a, [stim_start, stim_end],
        color=:gray, linestyle=:dash, alpha=0.5)

save("examples/plots/desensitization_off_bipolar_iGluR.png", fig3)
