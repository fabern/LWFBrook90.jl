# test script to set up Pfynwald model
include("LWFBrook90.jl")

input_prefix="pfynwald";
input_path="./examples/PFY2024-capacitance/";

model = LWFBrook90.loadSPAC(input_path, input_prefix; simulate_isotopes = false);

simulation = LWFBrook90.setup(model);

LWFBrook90.simulate!(simulation)

# extract output
sol = simulation.ODESolution;
timesteps = unique(round.(sol.t));
dates = LWFBrook90.RelativeDaysFloat2DateTime.(timesteps, simulation.parametrizedSPAC.reference_date);

PLSTOR = [sol(t).PLSTOR.mm for t in timesteps];
PLPSI = [sol(t).PLHYD.ψ for t in timesteps];
PLPSI_pd = [sol(t).accum.cum_pd_plpsi for t in timesteps];
PLPSI_md = [sol(t).accum.cum_md_plpsi for t in timesteps];
PLFL = [sol(t).accum.cum_d_plfl for t in timesteps];
RWU = [sol(t).RWU.mmday for t in timesteps];
TRANI = reduce(hcat, [sol(t).TRANI.mmday for t in timesteps]);
PLRFI = reduce(hcat, [sol(t).PLRFI.mmday for t in timesteps]);

soil = LWFBrook90.get_soil_(:ψ, simulation);
flux = LWFBrook90.get_fluxes(simulation);
flux[153:158, :]

#=
using Dates, DataFrames, Plots; gr()
plot(PLPSI_pd, label="Pre-dawn ψ")
plot(dates[153:158], PLPSI_pd[153:158], label="Pre-dawn ψ")
plot(122:303, PLPSI_pd[122:303], label="Pre-dawn ψ")
plot(dates, RWU, label="RWU")
plot(dates, PLFL, label="PLFL")
plot(dates, PLFL./(RWU.+PLFL), label="PLFL/RWU")

psi_comp = DataFrame(t=dates, PLPSI=PLPSI_pd)
psi_comp = [psi_comp; DataFrame(t=dates.+Hour(12), PLPSI=PLPSI_md)]
sort!(psi_comp, :t)
plot(psi_comp.t[305:315], psi_comp.PLPSI[305:315])

=#