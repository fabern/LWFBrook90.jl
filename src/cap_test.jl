# test script to set up Pfynwald model
include("LWFBrook90.jl")

input_prefix="pfynwald";
input_path="./examples/PFY2024-capacitance/";

model = LWFBrook90.loadSPAC(input_path, input_prefix; simulate_isotopes = true);

simulation = LWFBrook90.setup(model; requested_tspan = (152,157));

LWFBrook90.simulate!(simulation)

# extract output
sol = simulation.ODESolution;
timesteps = unique(round.(sol.t));

PLSTOR = [sol(t).PLSTOR.mm for t in timesteps];
PLPSI = [sol(t).PLHYD.ψ for t in timesteps];
PLPSI_pd = [sol(t).accum.cum_pd_plpsi for t in timesteps];
PLPSI_md = [sol(t).accum.cum_md_plpsi for t in timesteps];
PLFL = [sol(t).accum.cum_d_plfl for t in timesteps];
RWU = [sol(t).RWU.mmday for t in timesteps];
TRANI = [sol(t).TRANI.mmday for t in timesteps];

soil = LWFBrook90.get_soil_(:ψ, simulation);
flux = LWFBrook90.get_fluxes(simulation);

#=
using DataFrames, Plots; gr()
plot(PLPSI_pd)

psi_comp = DataFrame(t=timesteps, PLPSI=PLPSI_pd)
psi_comp = [psi_comp; DataFrame(t=timesteps.+0.5, PLPSI=PLPSI_md)]
sort!(psi_comp, :t)
plot(psi_comp.t, psi_comp.PLPSI)

=#