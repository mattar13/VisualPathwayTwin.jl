fig = Figure(size = (1200, 600))
ax1a = Axis(fig[1, 1], title = "Rho"); hidespines!(ax1a)
ax1b = Axis(fig[1, 2], title = "Tr"); hidespines!(ax1b)
ax2a = Axis(fig[2, 1], title = "PDE"); hidespines!(ax2a)
ax2b = Axis(fig[2, 2], title = "cGMP"); hidespines!(ax2b)
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
ax1 = Axis(fig2[1:2, 1], title = "V_t"); hidespines!(ax1)
ax2a = Axis(fig2[1, 2], title = "iPHOTO"); hidespines!(ax2a)
ax2b = Axis(fig2[1, 3], title = "iLEAK"); hidespines!(ax2b)
ax2c = Axis(fig2[1, 4], title = "iH"); hidespines!(ax2c)
ax2d = Axis(fig2[1, 5], title = "iKV"); hidespines!(ax2d)
ax3a = Axis(fig2[2, 2], title = "iCa"); hidespines!(ax2a)
ax3b = Axis(fig2[2, 3], title = "iKCa"); hidespines!(ax2b)
ax3c = Axis(fig2[2, 4], title = "iCl"); hidespines!(ax2c)

ax4a = Axis(fig2[3, 2], title = "iEX"); hidespines!(ax4a)
ax4b = Axis(fig2[3, 3], title = "iEX2"); hidespines!(ax4b)
ax4c = Axis(fig2[3, 4], title = "CA currents"); hidespines!(ax4c)
ax4d = Axis(fig2[3, 5], title = "Ca storage"); hidespines!(ax4d)


import DigitalTwin: hCa, mCl, mKCas, J∞
for (i, sol) in enumerate(data_series)
    J_t = map(t -> sol(t)[6], t_rng)
    V = map(t -> sol(t)[5], t_rng) #This is the voltage equation, but we are not using it in this model
    Gt = map(t -> sol(t)[4], t_rng)
    _Ca_s = map(t -> sol(t)[16], t_rng)
    O1 = map(t -> sol(t)[7], t_rng)
    O2 = map(t -> sol(t)[8], t_rng)
    O3 = map(t -> sol(t)[9], t_rng)
    mKV = map(t -> sol(t)[12], t_rng)
    hKV = map(t -> sol(t)[13], t_rng)
    mCa = map(t -> sol(t)[14], t_rng)
    mKCa = map(t -> sol(t)[15], t_rng)
    #Open parameters

    iPHOTO = @. -iDARK * J∞(Gt, kg)#* (1.0 - exp((V - 8.5) / 17.0))
    iLEAK = @. gLEAK*(V+eLEAK) #Leak
    iH =  @. gH*(O1 + O2 + O3)*(V + eH) #Ih Current
    iKV =  @. gKV*mKV^3+hKV*(V + eK)
    iCa =  @. gCa*mCa^4*hCa(V)*(V + eCa) #Ca current #We should add the log
    iKCa =  @. gKCa * mKCa^2 * mKCas(_Ca_s) * (V + eK) #KCa current
    iCl = @. gCl * mCl(_Ca_s) * (V + eCl) #Cl current
    iEX = @. J_ex * exp(-(V - 14) / 70) * (_Ca_s - Cae) / (_Ca_s - Cae + K_ex)
    iEX2 = @. J_ex2 * (_Ca_s - Cae) / (_Ca_s - Cae + K_ex2)

    lines!(ax2a, t_rng, iPHOTO,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))

    lines!(ax2b, t_rng, iLEAK,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))

    lines!(ax2c, t_rng, iH,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))

    lines!(ax2d, t_rng, iKV,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))

    lines!(ax3a, t_rng, iCa,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))

    lines!(ax3b, t_rng, iKCa,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))

    lines!(ax3c, t_rng, iCl,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))

    lines!(ax4a, t_rng, iEX,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))

    lines!(ax4b, t_rng, iEX2,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))

    lines!(ax4c, t_rng, J_t,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))

    V_t = map(t -> sol(t)[5], t_rng)
    lines!(ax1, t_rng, V_t,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))
end
fig2

#%%
 
(2*F*V1)

_Ca_0
gKCa
gCl
eCl
F
DCa
S1
DELTA
V1
S1 / (DELTA * V1)

#%%

fig3 = Figure(size = (1200, 600))
ax1 = Axis(fig3[1, 1], title = "Ca_s"); hidespines!(ax1)
ax2 = Axis(fig3[1, 2], title = "Ca_f"); hidespines!(ax2)
ax3 = Axis(fig3[2, 1], title = "CaB_ls"); hidespines!(ax3)
ax4 = Axis(fig3[3, 1], title = "CaB_hs"); hidespines!(ax4)
ax5 = Axis(fig3[2, 2], title = "CaB_lf"); hidespines!(ax5)
ax6 = Axis(fig3[3, 2], title = "CaB_hf"); hidespines!(ax6)
for (i, sol) in enumerate(data_series)

    CaS_t = map(t -> sol(t)[16], t_rng)
    CaF_t = map(t -> sol(t)[17], t_rng)
    CaB_ls_t = map(t -> sol(t)[18], t_rng)
    CaB_hs_t = map(t -> sol(t)[19], t_rng)
    CaB_lf_t = map(t -> sol(t)[20], t_rng)
    CaB_hf_t = map(t -> sol(t)[21], t_rng)

    lines!(ax1, t_rng, CaS_t,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))
    lines!(ax2, t_rng, CaF_t,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))
    lines!(ax3, t_rng, CaB_ls_t,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))
    lines!(ax4, t_rng, CaB_hs_t,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))
    lines!(ax5, t_rng, CaB_lf_t,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))
    lines!(ax6, t_rng, CaB_hf_t,     
        color = round(log10(photon_range[i])), colormap = :viridis, 
        colorrange = (0.0, 3.0))
end
fig3


#%% Make a second figure to show the gating of J
solE = data_series[end]
G_t = map(t -> solE(t)[4], t_rng)
A_t = map(t -> solE(t)[7], t_rng)
g_rng = LinRange(0.0, 2.0, 1000)
a_rng = LinRange(-1.0, 0.0, 1000)
Jinf_rng = map(g -> DigitalTwin.J∞(g, 20.0), g_rng)
Hinf_rng = map(a -> DigitalTwin.H_inf(a, 0.000898898, 0.020486497), a_rng)

fig2 = Figure(size = (1200, 600))
ax1 = Axis(fig2[1, 1], title = "G"); hidespines!(ax1)
ax2 = Axis(fig2[2, 1], title = "H"); hidespines!(ax2)
ax3 = Axis(fig2[1, 2], title = "G"); hidespines!(ax3)
ax4 = Axis(fig2[2, 2], title = "A"); hidespines!(ax4)

lines!(ax1, t_rng, G_t, color = :red, label = "J")
lines!(ax2, t_rng, -30.0*map(g -> DigitalTwin.J∞(g, 0.1), G_t), color = :red, label = "H")

lines!(ax3, g_rng, Jinf_rng, color = :red, label = "H")
lines!(ax4, a_rng, Hinf_rng, color = :red, label = "H")

fig2