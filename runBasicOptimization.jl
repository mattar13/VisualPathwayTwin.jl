using Revise
using Pkg; Pkg.activate(".")
using ElectroPhysiology, PhysiologyModeling
#import ElectroPhysiology.addStimulus!
using DiffEqParamEstim, Optimization
using Statistics
using Optim #Add to main package
#Use these packages for plotting and saving
using Pkg;Pkg.activate("test")
using GLMakie, PhysiologyPlotting
using DataFrames, CSV

#%% Input the working data
include("OpenData.jl")
fn = raw"E:\Data\ERG\Retinoschisis\2021_05_27_RS1KO-30\Mouse1_Adult_RS1KO\NoDrugs\Rods\nd1_1p_0000.abf"
fn = raw"E:\Data\ERG\Melanopsin Data\2022_04_21_MelCreAdult\Mouse2_Adult_MelCre\NoDrugs\Rods\nd1_1p_0000.abf"
fn = raw"E:\Data\ERG\Melanopsin Data\2022_04_21_MelCreAdult\Mouse2_Adult_MelCre\BaCl_LAP4\Rods\nd1_1p_0000.abf"

dataERG, expERG = openData(fn)
dataERG.t
stim_start = getStimulusStartTime(dataERG)[1]*1000
stim_end = getStimulusEndTime(dataERG)[1]*1000
Stim(t) = stim_start <= t <= stim_end ? 400.0 : 0.0
println("Stimulus runs from $(stim_start) to $(stim_end)")

getchannel(dataERG, 1).data_array[1,:,1]

#%% Set up the optimization problem
include("Models.jl")
param_df = CSV.read(raw"E:\KozLearn\Standards\starting_params.csv", DataFrame)
keys = param_df.Key
p0 = param_df.Value
lower_bounds = param_df.LowerBounds
upper_bounds = param_df.UpperBounds

#%% Maybe we can try to setup black box optim as a initial guess
a_p0 = p0[[1:9..., 19:23...]]

#using BlackBoxOptim
# --- Optimization settings---
opt = Optim.Options(
    #g_abstol = 1e-8, 
    #g_reltol = 1e-8,
    iterations = 100,
    show_trace = true, 
    show_every = 10
)

alg = Fminbox(BFGS())
# Use Optim.jl to minimize the loss function with respect to parameters.
result = optimize(y -> loss_static(dataERG, y), 
    lower_bounds, upper_bounds, p0, 
    alg, opt, autodiff = :forward
)

#Get the results
opt_params = Optim.minimizer(result)
# Extract the optimized parameters and print them.
println("Optimized parameters: ", opt_params, "with loss of: ", loss_static(dataERG, opt_params))

#%%
#opt_params = p0
#opt_params = [0.3932316550412567, 1.0260082770989687, 0.0003916972931226907, 800.0006939405749, 0.8939308786851997, 10.007262117944938, 100.04652651427836, 0.02854849213840683, 0.09999999972096212, 59.99999999999999, 3.328106018630575e-6, 3.5094030906847675, 9.559471340694554e-6, 0.1974602085556033, 0.5044543581561698, 0.5325954703755535, 1.0150865919213672, 0.2087574134432508, 14.866287408516474, 9.697855871116372, 38.8862854460252, 240.86751760640078, 44.68145218102656, 134.30804769035677, 519.2998392601478, 31816.756356838676, 199.9996841782422]

#%% Process the plotting
stim_t = map(x -> Stim(x), dataERG.t)
sol, ERG_t = simulate_awave(dataERG.t, opt_params[[1:9..., 19:23...]], opt_params)
sol_t = sol.t
#Plot components
a_t = map(t -> sol(t)[6], sol_t)
b_t = map(t -> sol(t)[7], sol_t)
m_t = map(t -> sol(t)[8], sol_t)
c_t = map(t -> sol(t)[9], sol_t)
o_t = map(t -> sol(t)[13], sol_t)
#ERG_t = map(t->sol(t)[14], sol_t)

# Plot the results
fig = Figure(size = (1000, 600))
ax1 = Axis(fig[1, 1], title = "ERG data collected"); hidedecorations!(ax1); hidespines!(ax1)
experimentplot!(ax1, dataERG, channel = 3)
vlines!(ax1, [stim_end], color = :green, label = "Stimulus", alpha = 0.5)
ylims!(ax1, (-1.0, 0.20))
#Add all of the components together

sb_x = 2500
sb_y = 0.2
sb_dx = 250
sb_dy = 0.25
lines!(ax1, [sb_x, sb_x+sb_dx], [sb_y, sb_y], color = :black, linewidth = 3)  # horizontal scale bar
lines!(ax1, [sb_x, sb_x], [sb_y, sb_y+sb_dy], color = :black, linewidth = 3)  # vertical scale bar
text!(ax1, sb_x-(sb_dx/3), sb_y + (sb_dy/2), text = "$(ceil(Int64, sb_dy*1000))μV", align = (:center, :center), rotation = π/2, color = :black, fontsize = 12.0, font = :bold)
text!(ax1, sb_x+(sb_dx/2), sb_y - (sb_dy/3), text = "$(ceil(Int64, sb_dx)) ms", align = (:center, :center), color = :black, fontsize = 12.0, font = :bold)

ax2 = Axis(fig[1, 2], title = "Modelled Retina"); hidedecorations!(ax2); hidespines!(ax2)
experimentplot!(ax2, dataERG, channel = 3)
lines!(ax2, sol_t, ERG_t, color = :red, label = "Simulated ERG")
vlines!(ax2, [stim_end], color = :green, label = "Stimulus", alpha = 0.5)

sb_x = 2500
sb_y = 0.2
sb_dx = 250
sb_dy = 0.25
lines!(ax1, [sb_x, sb_x+sb_dx], [sb_y, sb_y], color = :black, linewidth = 3)  # horizontal scale bar
lines!(ax1, [sb_x, sb_x], [sb_y, sb_y+sb_dy], color = :black, linewidth = 3)  # vertical scale bar
text!(ax1, sb_x-(sb_dx/3), sb_y + (sb_dy/2), text = "$(ceil(Int64, sb_dy*1000))μV", align = (:center, :center), rotation = π/2, color = :black, fontsize = 12.0, font = :bold)
text!(ax1, sb_x+(sb_dx/2), sb_y - (sb_dy/3), text = "$(ceil(Int64, sb_dx)) ms", align = (:center, :center), color = :black, fontsize = 12.0, font = :bold)

ax3 = Axis(fig[2,2]) ; hidedecorations!(ax3); hidespines!(ax3)
lines!(ax3, [NaN], [NaN], color = :red, label = "Simulated ERG")
lines!(ax3, [NaN], [NaN], color = :black, label = "Collected ERG")
vlines!(ax3, [stim_end], linewidth = 2.0, color = :green, label = "Stimulus (1ms bright light)", alpha = 0.5)
lines!(ax3, sol_t, a_t, linewidth = 2.0, color = :orange, label = "a-wave")
lines!(ax3, sol_t, b_t, linewidth = 2.0, color = :blue, label = "b-wave")
lines!(ax3, sol_t, m_t, linewidth = 2.0, color = :darkorchid, label = "Müller component")
lines!(ax3, sol_t, c_t, linewidth = 2.0, color = :darkcyan, label = "c-wave")
lines!(ax3, sol_t, o_t, color = :darkred, label = "OPs")
linkxaxes!(ax1, ax2, ax3)
linkyaxes!(ax1, ax2, ax3)
Legend(fig[2,1], ax3, tellwidth = false, halign = :right)

fig
#%%
 # Plot the results as a vector
y = [11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 11, 10, 9, 8, 7, 6]
x = fill(1, length(y))
x[12:end] .= 3

ax4 = Axis(fig[1:2, 3], ); hidespines!(ax4); hidedecorations!(ax4)
labels = keys
pars = opt_params
marker_colors = [
    :orange, 
    :orange, 
    :blue, 
    :darkorchid, 
    :darkcyan, 
    :darkorchid, 
    :darkred, 
    :darkred, 
    :darkred,
    :darkred, 
    :darkred, 
    :orange, 
    :orange, 
    :blue, 
    :darkorchid, 
    :darkcyan,
    :darkred

]
scatter!(ax4, x, y;
    #transparency = true, 
    markersize = 52,
    marker = :circle,
    alpha = 0.6,
    color = log10.(pars),
    colormap = :viridis,
    strokecolor = marker_colors,        # dark outline
    strokewidth = 3)

for (i, lab) in enumerate(labels)
    text!(ax4, lab, position = (x[i]-0.75, y[i]), align = (:center, :center), fontsize = 20, color = :black)
    if lab == "tC"
        num = "$(round(pars[i]/10000, digits = 2))e5"
    elseif lab == "tM"
        num = "$(round(pars[i]/100, digits = 2))e3"
    else
        num = "$(round(pars[i], digits = 2))"
    end
    text!(ax4, num, position = (x[i], y[i]), align = (:center, :center), fontsize = 12, color = :black)

end
colsize!(fig.layout, 3, Relative(1/4))
xlims!(ax4, -0.0, 4.0)

fig

#%% Save the plot
save(raw"E:\KozLearn\Standards\no_drugs_rs1_trace.png", fig)

#%% Save the data if you like it
# Create a DataFrame to hold the parameter index and value
df = DataFrame(Key = keys, Value = opt_params, LowerBounds = lower_bounds, UpperBounds = upper_bounds)

# Write the DataFrame to a CSV file
CSV.write(raw"E:\KozLearn\Standards\no_drugs_params.csv", df)
#CSV.write(raw"E:\KozLearn\Standards\standard_params.csv", df)