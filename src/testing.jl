using Revise
using DigitialTwin

using GLMakie, PhysiologyPlotting


using DataFrames, CSV

stim_start = 10.0
stim_end = 11.0
photons = 400.0

#make the model

params = [3.5]
u = zeros(8)
tspan = (0.0, 100.0)
model = ODEProblem(erg_ode!, u, tspan, params)
sol = solve(model, Tsit5(), saveat = 0.1)

