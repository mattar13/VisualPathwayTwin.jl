fig = Figure(size = (1200, 600))
ax1a = Axis(fig[1, 1], title = "Rho"); hidespines!(ax1a)
ax1b = Axis(fig[1, 2], title = "Tr"); hidespines!(ax1b)
ax2a = Axis(fig[2, 1], title = "PDE"); hidespines!(ax1c)
ax2b = Axis(fig[2, 2], title = "cGMP"); hidespines!(ax1d)
ax3 = Axis(fig[1:2, 3], title = "iPHOTO"); hidespines!(ax3)

for (i, sol) in enumerate(data_series)
    R_t = map(t -> sol(t)[1], t_rng)
    T_t = map(t -> sol(t)[2], t_rng)
    P_t = map(t -> sol(t)[3], t_rng)
    G_t = map(t -> sol(t)[4], t_rng)
    J_t = map(t -> sol(t)[6], t_rng)

    vlines!(ax1a, [stim_start, stim_end], (-0.5, 1.0), color = :black, alpha = 0.2)
    lines!(ax1a, t_rng, R_t,         
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0),
        label = "Tr")

    vlines!(ax1b, [stim_start, stim_end], (-0.5, 1.0), color = :black, alpha = 0.2)
    lines!(ax1b, t_rng, T_t,         
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0),
        label = "Tr")

    vlines!(ax2a, [stim_start, stim_end], (-0.5, 1.0), color = :black, alpha = 0.2)
    lines!(ax2a, t_rng, P_t,         
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0),
        label = "Tr")

    vlines!(ax2b, [stim_start, stim_end], (-0.5, 1.0), color = :black, alpha = 0.2)
    lines!(ax2b, t_rng, G_t,         
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0),
        label = "Tr")

    vlines!(ax3, [stim_start, stim_end], (-0.5, 1.0), color = :black, alpha = 0.2)
    lines!(ax3, t_rng, J_t,         
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0), label = "A-wave")
end
fig

#%% Make a second figure to show the gating of J
fig2 = Figure(size = (1200, 600))
ax1a = Axis(fig2[1, 1], title = "Rho"); hidespines!(ax1a)
ax1b = Axis(fig2[1, 2], title = "Tr"); hidespines!(ax1b)
ax2a = Axis(fig2[2, 1], title = "PDE"); hidespines!(ax1c)
ax2b = Axis(fig2[2, 2], title = "cGMP"); hidespines!(ax1d)
ax3 = Axis(fig2[1:2, 3], title = "iPHOTO"); hidespines!(ax3)

for (i, sol) in enumerate(data_series)
    J_t = map(t -> sol(t)[6], t_rng)







    V_t = map(t -> sol(t)[5], t_rng)
    lines!(ax3, t_rng, V_t,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0), label = "V")
end
fig2

#%%
HC1_t = map(t -> sol(t)[7], t_rng)
HC2_t = map(t -> sol(t)[8], t_rng)
HO1_t = map(t -> sol(t)[9], t_rng)
HO2_t = map(t -> sol(t)[10], t_rng)
HO3_t = map(t -> sol(t)[11], t_rng)

CaS_t = map(t -> sol(t)[16], t_rng)
CaF_t = map(t -> sol(t)[17], t_rng)
CaB_ls_t = map(t -> sol(t)[18], t_rng)
CaB_hs_t = map(t -> sol(t)[19], t_rng)
CaB_lf_t = map(t -> sol(t)[20], t_rng)
CaB_hf_t = map(t -> sol(t)[21], t_rng)
