# test script to set up Pfynwald model
include("LWFBrook90.jl")

input_prefix="pfynwald";
input_path="./examples/PFY2024-capacitance/";

model = LWFBrook90.loadSPAC(input_path, input_prefix; simulate_isotopes = false,
    root_distribution = (beta = 0.97091, z_rootMax_m=-1.35462));

simulation = LWFBrook90.setup(model, requested_tspan = (0, 3));

LWFBrook90.simulate!(simulation)

# extract output
sol = simulation.ODESolution;
timesteps = unique(round.(sol.t));
dates = LWFBrook90.RelativeDaysFloat2DateTime.(timesteps, simulation.parametrizedSPAC.reference_date);

SWATI = [sol(t).SWATI.mm for t in timesteps];
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

# post-processing of new saved_all

out = simulation.saved_all;
t = out.t;

RWU_fix = [out.saveval[t].RWU for t in 1:length(t)];
prec_fix = [out.saveval[t].accum.cum_d_prec for t in 1:length(t)];

# compare accumulated meteo with forcing
meteo = LWFBrook90.get_forcing(simulation);
all((meteo.prec_mmDay[1:(length(t)-1)] .- prec_fix[1:(length(t)-1)]) .< 1e-3);

using DataFrames

function saved_values_to_dataframe(saved)
    nt = length(saved.t)
    all_keys = fieldnames(typeof(saved.saveval[1]))  # keys of named tuple

    data_cols = Dict{Symbol, Any}()
    data_cols[:t] = saved.t

    for k in all_keys
        # Extract all values for this key over time
        vals = [getfield(saved.saveval[i], k) for i in 1:nt]

        if (vals[1] isa AbstractArray)
            nsub = length(vals[1])
            colnames = propertynames(vals[1])
            for j in 1:nsub
                if length(colnames) > 1
                    colname = colnames[j]
                else
                    colname = Symbol("$(k)_$j")
                end
                data_cols[colname] = [v[j] for v in vals]
            end
        else
            data_cols[Symbol(k)] = vals
        end
    end

    return DataFrame(data_cols)
end

flux_fix = saved_values_to_dataframe(simulation.saved_all)