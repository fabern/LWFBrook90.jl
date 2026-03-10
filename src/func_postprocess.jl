#######
####### Public facing API: getting simulation results ======================================
#######

@doc raw"""
    get_states(simulation::DiscretizedSPAC; days_to_read_out_d = nothing)

Returns a DataFrame with amounts (and isotopoic compositions) of the inputs and state variables:
PREC, GWAT, INTS, INTR, SNOW, (amount only: CC, SNOWLQ), (signature only: XYLEM). RWU flux can be obtained with `get_fluxes()`.
By default, the values are returned for each simulation day.
The user can optionally define timestep using the input argument `days_to_read_out_d` as:
- `:integrator_step` to have it each simulation time step
- a numeric vector, e.g. `days_to_read_out_d = 1:1.0:100` for specific days since `simulation.parametrizedSPAC.reference_date`

Subsets of scalar state variables or soil states can be mad afterwards with:

    simout_states = get_states(simulation)
    simout_states[:, Not(r"[0-9]mm$")]    # select only scalar states (by de-selecting columns containing depth information e.g. "_150mm")
    simout_states[:, r"(dates)|([0-9]mm)"] # select only vector states (by selecting columns containing depth information e.g. "_150mm")

Function returns a DataFrame:

    366×108 DataFrame
     Row │ dates                GWAT_mm  INTS_mm   INTR_mm       SNOW_mm  CC_MJm2     SNOWLQ_mm  SWAT_mm  GWAT_d18O  INTS_d18O  INTR_d18O  SNOW_d18O  XYLEM_d18O  ⋯ δ18O_permil_50mm  δ18O_permil_100mm  δ18O_permil_1 ⋯ ψ_kPa_1100mm
         │ DateTime             Float64  Float64   Float64       Float64  Float64     Float64    Float64  Float64    Float64    Float64    Float64    Float64     ⋯ Float64           Float64            Float64       ⋯ Float 64
    ─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
       1 │ 2021-01-01T00:00:00      1.0  0.0        0.0           0.0     0.0          0.0       90.2641  -13.0        -13.0      -13.0      -13.0     -10.1111  ⋯         -10.1111           -10.1111           -10. ⋯ -7
       2 │ 2021-01-02T00:00:00      1.0  0.0        0.0           0.0     0.0          0.0       90.4532  -12.5973     -15.04     -15.04     -15.04    -10.1111            -10.1111           -10.1111           -10. ⋯ -6.43
"""
function get_states(simulation::DiscretizedSPAC; days_to_read_out_d = nothing) # returns scalar states GWAT, INTS, INTR, SNOW, CC (amount only), SNOWLQ (amount only), XYLEM (signature only)
    @assert !isnothing(simulation.ODESolution) "Solution was not yet computed. Please simulate!(simulation)"
    timepoints =
        isnothing(days_to_read_out_d)           ? unique(round.(simulation.ODESolution.t)) : # case nothing: return daily by default
        days_to_read_out_d == :integrator_step  ? simulation.ODESolution.t :                 # if requested return for each integration step stored in ODESolution
        days_to_read_out_d[1] isa Number        ? days_to_read_out_d :                       # else assume we got a vector of days
        error("Unknown `days_to_read_out_d` provided.")

    # i) scalar states (+PREC forcing)
    simulate_isotopes = simulation.parametrizedSPAC.solver_options.simulate_isotopes

    states = hcat(
        LWFBrook90.intern___get_scalars(
            [:GWAT, :INTS, :INTR, :SNOW, :CC,   :SNOWLQ, :SWATI],
            [:mm,   :mm,   :mm,   :mm,   :MJm2, :mm,     :mm],
            simulation, timepoints),
        !simulate_isotopes ? DataFrame() : LWFBrook90.intern___get_scalars(
            [:GWAT, :INTS, :INTR, :SNOW, :XYLEM], # , :RWU
            [:d18O, :d18O, :d18O, :d18O, :d18O],  # , :d18O
            simulation, timepoints)[:,Not(:time)],
        !simulate_isotopes ? DataFrame() : LWFBrook90.intern___get_scalars(
            [:GWAT, :INTS, :INTR, :SNOW, :XYLEM], # , :RWU
            [:d2H,  :d2H,  :d2H,  :d2H,  :d2H],   # , :d2H
            simulation, timepoints)[:,Not(:time)]) # alternatively innerjoin on = [:time])
    rename!(states, :SWATI_mm => :SWAT_mm)
    states.SWAT_mm = sum.(states.SWAT_mm)
    # scalar states:
    # simulation.ODESolution.u[1].GWAT
    # simulation.ODESolution.u[1].INTS
    # simulation.ODESolution.u[1].INTR
    # simulation.ODESolution.u[1].SNOW
    # simulation.ODESolution.u[1].XYLEM
    # simulation.ODESolution.u[1].CC.MJm2
    # simulation.ODESolution.u[1].SNOWLQ.mm

    # ii) vector states
    # names(get_soil_([:θ, :ψ, :δ18O, :δ2H, :W, :SWATI, :K], simulation))
    vect_states = get_soil_([:θ, :ψ, :d18O, :d2H], simulation; days_to_read_out_d = timepoints)

    # simulation.ODESolution.u[1].SWATI.mm
    # simulation.ODESolution.u[1].SWATI.d2H
    # simulation.ODESolution.u[1].SWATI.d18O
    # simulation.ODESolution.u[1].aux.θ
    # simulation.ODESolution.u[1].aux.ψ
    # simulation.ODESolution.u[1].aux.K

    # Define dates as DateTime from time:
    t_ref = simulation.ODESolution.prob.p.REFERENCE_DATE
    # dates = RelativeDaysFloat2DateTime.(days_to_read_out_d, t_ref)
    states.time      = RelativeDaysFloat2DateTime.(states.time,      t_ref); rename!(states,      :time => :dates)
    vect_states.time = RelativeDaysFloat2DateTime.(vect_states.time, t_ref); rename!(vect_states, :time => :dates)

    return innerjoin(states, vect_states, on = [:dates])
end

@doc raw"""
    get_fluxes(simulation::DiscretizedSPAC; days_to_read_out_d = nothing)

Returns a DataFrame with amounts (and isotopoic compositions) of the ecosystem water fluxes*.

By default, the values are returned for each simulation day.
The user can optionally define timestep using the input argument `days_to_read_out_d` as:
- a numeric vector, e.g. `days_to_read_out_d = 1:1.0:100` for specific days since `simulation.parametrizedSPAC.reference_date`

Subsets of scalar state variables or soil states can be created afterwards with:

    simout_fluxes = get_fluxes(simulation)
    simout_fluxes[:, Not(r"[0-9]mm$")]    # select only scalar fluxes (by de-selecting columns containing depth information e.g. "_150mm")
    simout_fluxes[:, r"(dates)|([0-9]mm)"] # select only vector fluxes (by selecting columns containing depth information e.g. "_150mm")

Naming of fluxes is the same as in the illustration in the documentation. (TODO, copy paste figure and table from article into documentation.)
Not shown in illustration are PINT, PTRAN, PSLVP, which correspond to the potential
interception rate (i.e. total amount of evaporated water from a continuously wet canopy),
potential transpiration rate, and potential ground evaporation rate (of soil water or snow).

Function returns a DataFrame:

    366×50 DataFrame
     Row │ dates                cum_d_prec  cum_d_sfal  cum_d_sthr  cum_d_sint  cum_d_rfal  cum_d_rint  cum_d_rthr  cum_d_rsno  cum_d_rnet  cum_d_smlt   ⋯
         │ DateTime             Float64     Float64     Float64     Float64     Float64     Float64     Float64     Float64     Float64     Float64      ⋯
    ─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── ⋯
       1 │ 2021-01-01T00:00:00         0.0    0.0         0.0        0.0          0.0         0.0         0.0         0.0          0.0         0.0       ⋯
       2 │ 2021-01-02T00:00:00         0.0    0.0         0.0        0.0          0.0         0.0         0.0         0.0          0.0         0.0       ⋯

         ⋯ cum_d_irvp  cum_d_isvp  cum_d_snvp  cum_d_slvp  cum_d_tran  ⋯
         ⋯ Float64     Float64     Float64     Float64     Float64     ⋯
         ⋯ ─────────────────────────────────────────────────────────── ⋯
         ⋯   0.0        0.0         0.0        0.0         0.0         ⋯
         ⋯   0.0        0.0         0.0        0.0448292   0.0566162   ⋯

         ⋯ cum_d_pint  cum_d_ptran  cum_d_pslvp ⋯ srfl     slfl      byfl     dsfl      gwfl      vrfln     ⋯
         ⋯ Float64     Float64      Float64     ⋯ Float64  Float64   Float64  Float64   Float64   Float64   ⋯
         ⋯ ──────────────────────────────────────────────────────────────────────────────────────────────── ⋯
         ⋯ 0.0       0.0            0.0         ⋯   0.0   0.0          0.0      0.0   0.0       0.0         ⋯
         ⋯ 2.64652   0.0566162      0.364727    ⋯   0.0   0.0          0.0      0.0   0.154122  0.154122    ⋯

         ⋯ flow      seep     evap       TRANI_mmday_50mm  TRANI_mmday_100mm  TRANI_mmday_150mm  TRANI_mmday_200mm  TRANI_mmday_250mm  ⋯ TRANI_mmday_1050mm TRANI_mmday_1100mm
         ⋯ Float64   Float64  Float64    Float64           Float64            Float64            Float64            Float64            ⋯ Float64            Float64
         ⋯ ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
         ⋯ 0.0           0.0  0.0             0.0                0.0                0.0                0.0                0.0          ⋯                0.0               0.0
         ⋯ 0.154122      0.0  0.101445        0.0434988          0.00995972         0.00240976         0.000574825        0.00013463   ⋯                0.0               0.0
"""
function get_fluxes(simulation::DiscretizedSPAC; days_to_read_out_d = nothing) # returns fluxes (external, internal, errorterms)
    @assert !isnothing(simulation.ODESolution) "Solution was not yet computed. Please simulate!(simulation)"

    timepoints =
        isnothing(days_to_read_out_d)           ? unique(round.(simulation.ODESolution.t)) : # case nothing: return daily by default
        # days_to_read_out_d == :integrator_step  ? error("Cumulative fluxes should not be read out on subdaily intervals.") : # if requested return for each integration step stored in ODESolution
        days_to_read_out_d == :integrator_step  ? simulation.ODESolution.t :
        days_to_read_out_d[1] isa Number        ? days_to_read_out_d :                       # else assume we got a vector of days
        error("Unknown `days_to_read_out_d` provided.")
    if !isnothing(days_to_read_out_d) @warn("""
        Provided `days_to_read_out_d` to `get_fluxes()`.
        Be careful at what points you read out the average or daily cumulative fluxes to use them properly.
        It might be simpler to read out daily fluxes (by not specifiying `days_to_read_out_d`) and aggregate later with more control. """)
    end

    t_ref = simulation.ODESolution.prob.p.REFERENCE_DATE

    # 1) scalar fluxes
    ks = keys(simulation.ODESolution.u[1].accum)
    df_scalar_fluxes =
        DataFrame((:dates => RelativeDaysFloat2DateTime.(timepoints, t_ref)),
                (k => [simulation.ODESolution(t).accum[k] for t in timepoints] for k in ks)...)
    # DataFrame((:date => simulation.ODESolution_datetime ),
    #         (k => [ut.accum[k] for ut in simulation.ODESolution.u] for k in ks)...)

    # 2) vector fluxes
    # trani:
    df_vector_fluxes = get_soil_([:TRANI], simulation; days_to_read_out_d = timepoints)
    # TODO: # byfli, infli, dsfli are currently not done

    t_ref = simulation.ODESolution.prob.p.REFERENCE_DATE
    df_vector_fluxes.time = RelativeDaysFloat2DateTime.(df_vector_fluxes.time, t_ref)
    rename!(df_vector_fluxes, :time => :dates)

    # 3) scalar fluxes, signatures
    simulate_isotopes = simulation.parametrizedSPAC.solver_options.simulate_isotopes
    df_scalar_signatures = !simulate_isotopes ? DataFrame() : LWFBrook90.intern___get_scalars(
        [:RWU, :PREC, :RWU, :PREC, :RWU],
        [:mmday,:d18O, :d18O, :d2H, :d2H],
        simulation, timepoints)[:,Not(:time)]

    # 4) vector fluxes, signatures
    df_vector_signatures = DataFrame() # TODO: this is currently not stored anywhere when simulating (would need to append to state vector u0)

    # combine scalar and vector
    df_all_fluxes = innerjoin(df_scalar_fluxes, df_vector_fluxes, on = [:dates])
    df_all_fluxes = hcat(df_all_fluxes, df_scalar_signatures, df_vector_signatures)

    # reorder columns
    if simulate_isotopes
        select!(df_all_fluxes,
        :dates,
        # i) internal fluxes:
        :cum_d_prec,:cum_d_sfal,:cum_d_sthr,:cum_d_sint,:cum_d_irrig,
        :cum_d_rfal,:cum_d_rint,:cum_d_rthr,:cum_d_rsno,:cum_d_rnet,:cum_d_smlt,
        :cum_d_irvp,:cum_d_isvp,:cum_d_snvp,:cum_d_slvp,
        :cum_d_plfl,:cum_d_plrf,:cum_pd_plpsi,:cum_md_plpsi,
        :cum_d_tran, # TODO: :accum.cum_d_tran, # delete. we can use u[1].RWU.mmday # all([(simulation.ODESolution.u[idx].accum.cum_d_tran == simulation.ODESolution.u[idx].RWU.mmday) for idx in eachindex(simulation.ODESolution.u)])
        # :RWU_mmday, # use RWU instead of accum.cum_d_tran
        :RWU_d18O, :RWU_d2H,
        :PREC_d18O, :PREC_d2H,

        :cum_d_pint, :cum_d_ptran, :cum_d_pslvp, # not shown in Figure of P2
        :srfl,:slfl,:byfl,:dsfl,:gwfl,:vrfln,

        # ii) external fluxes:
        # :cum_d_prec
        :flow, :seep, :evap, # formerly cum_d_evap

        # # iii) for error computation
        # simulation.ODESolution.u[1].accum.StorageSWAT
        # simulation.ODESolution.u[1].accum.StorageWATER
        # simulation.ODESolution.u[1].accum.BALERD_SWAT
        # simulation.ODESolution.u[1].accum.BALERD_total

        # iv) internal fluxes per soil layer
        # # # TODO: add byfli # define this as accum.vec_byfl.mmday, .d18O, .d2H # similar to u.TRANI
        # # # TODO: add infli # define this as accum.vec_infl.mmday, .d18O, .d2H # similar to u.TRANI
        # # # TODO: add dsfli # define this as accum.vec_dsfl.mmday, .d18O, .d2H # similar to u.TRANI
        # # simulation.ODESolution.u[1].TRANI.mmday
        r"TRANI_"
        )
    else
        select!(df_all_fluxes,
        :dates,
        # i) internal fluxes:
        :cum_d_prec,:cum_d_sfal,:cum_d_sthr,:cum_d_sint,:cum_d_irrig,
        :cum_d_rfal,:cum_d_rint,:cum_d_rthr,:cum_d_rsno,:cum_d_rnet,:cum_d_smlt,
        :cum_d_irvp,:cum_d_isvp,:cum_d_snvp,:cum_d_slvp,
        :cum_d_plfl,:cum_d_plrf,:cum_pd_plpsi,:cum_md_plpsi,
        :cum_d_tran, # TODO: :accum.cum_d_tran, # delete. we can use u[1].RWU.mmday # all([(simulation.ODESolution.u[idx].accum.cum_d_tran == simulation.ODESolution.u[idx].RWU.mmday) for idx in eachindex(simulation.ODESolution.u)])
        # :RWU_mmday, # use RWU instead of accum.cum_d_tran

        :cum_d_pint, :cum_d_ptran, :cum_d_pslvp, # not shown in Figure of P2
        :srfl,:slfl,:byfl,:dsfl,:gwfl,:vrfln,

        # ii) external fluxes:
        # :cum_d_prec
        :flow, :seep, :evap, # formerly cum_d_evap

        # # iii) for error computation
        # simulation.ODESolution.u[1].accum.StorageSWAT
        # simulation.ODESolution.u[1].accum.StorageWATER
        # simulation.ODESolution.u[1].accum.BALERD_SWAT
        # simulation.ODESolution.u[1].accum.BALERD_total

        # iv) internal fluxes per soil layer
        # # # TODO: add byfli # define this as accum.vec_byfl.mmday, .d18O, .d2H # similar to u.TRANI
        # # # TODO: add infli # define this as accum.vec_infl.mmday, .d18O, .d2H # similar to u.TRANI
        # # # TODO: add dsfli # define this as accum.vec_dsfl.mmday, .d18O, .d2H # similar to u.TRANI
        # # simulation.ODESolution.u[1].TRANI.mmday
        r"TRANI_"
        )
    end

    return df_all_fluxes

    # Reorder columns
    #= return select(df_all_fluxes,
        :dates,
        # i) internal fluxes:
        :cum_d_prec,:cum_d_sfal,:cum_d_sthr,:cum_d_sint,:cum_d_irrig,
        :cum_d_rfal,:cum_d_rint,:cum_d_rthr,:cum_d_rsno,:cum_d_rnet,:cum_d_smlt,
        :cum_d_irvp,:cum_d_isvp,:cum_d_snvp,:cum_d_slvp,
        :cum_d_plfl,:cum_d_plrf,:cum_pd_plpsi,:cum_md_plpsi,
        :cum_d_tran, # TODO: :accum.cum_d_tran, # delete. we can use u[1].RWU.mmday # all([(simulation.ODESolution.u[idx].accum.cum_d_tran == simulation.ODESolution.u[idx].RWU.mmday) for idx in eachindex(simulation.ODESolution.u)])
        # :RWU_mmday, # use RWU instead of accum.cum_d_tran
        :RWU_d18O, :RWU_d2H,
        :PREC_d18O, :PREC_d2H,

        :cum_d_pint, :cum_d_ptran, :cum_d_pslvp, # not shown in Figure of P2
        :srfl,:slfl,:byfl,:dsfl,:gwfl,:vrfln,

        # ii) external fluxes:
        # :cum_d_prec
        :flow, :seep, :evap, # formerly cum_d_evap

        # # iii) for error computation
        # simulation.ODESolution.u[1].accum.StorageSWAT
        # simulation.ODESolution.u[1].accum.StorageWATER
        # simulation.ODESolution.u[1].accum.BALERD_SWAT
        # simulation.ODESolution.u[1].accum.BALERD_total

        # iv) internal fluxes per soil layer
        # # # TODO: add byfli # define this as accum.vec_byfl.mmday, .d18O, .d2H # similar to u.TRANI
        # # # TODO: add infli # define this as accum.vec_infl.mmday, .d18O, .d2H # similar to u.TRANI
        # # # TODO: add dsfli # define this as accum.vec_dsfl.mmday, .d18O, .d2H # similar to u.TRANI
        # # simulation.ODESolution.u[1].TRANI.mmday
        r"TRANI_"
        ) =#
                            # df_partitioning_daily = @chain df_fluxes begin
                            #         @rtransform begin
                            #             :ETa           = :evap # is actual evapotranspiration, i.e. sum of IRVP + ISVP + SNVP + SLVP + sum(aux_du_TRANI)
                            #             :Esoil         = :cum_d_slvp
                            #             :Esnow         = :cum_d_snvp
                            #             :Einterception = :cum_d_irvp + :cum_d_isvp
                            #             :Ta            = :cum_d_tran
                            #             :Precip        = :cum_d_prec
                            #             :Td            = :cum_d_ptran - :cum_d_tran
                            #             :D             = - :vrfln
                            #             # :R1          = -(:flow - :vrfln)    # flow = srfl+byfl+dsfli+gwfl, gwfl, vrfln
                            #             :R             = -(:srfl + :byfl + :dsfl) # This is more correctly not accounting for gwfl, thereby excluding state variable GWAT from the balance
                            #             :Swat          = :StorageSWAT
                            #         end
                            #         @select(:date, :year, :month, :ETa,:Esoil,:Esnow,:Einterception,:Ta,:Precip,:Td,:D,:R,:Swat)
                            #     end
end

@doc raw"""
    get_forcing(simulation::DiscretizedSPAC)

Returns a DataFrame with amounts (and isotopoic compositions) of the input forcing.
Note that meteorologic forcing (amounts and isotopic signatures) is same as in input csvs,
vegetation evolution (lai, height, sai) have been transformed to absolute units.

     Row │ dates                globrad_MJDayM2  tmax_degC  tmin_degC  vappres_kPa  windspeed_ms  prec_mmDay  precdelta18O_permil  precdelta2H_permil  densef_percent  height_m  lai_m2m2  sai_m2m2
         │ DateTime             Float64          Float64    Float64    Float64      Float64       Float64     Float64              Float64             Float64         Float64   Float64   Float64
    ─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
       1 │ 2021-01-01T00:00:00             5.53        0.9      -10.1         0.26           1.5         0.0               -15.04             -111.96           100.0      25.0       2.1       1.0
       2 │ 2021-01-02T00:00:00             5.1        -0.6       -9.3         0.3            1.6         0.0               -15.04             -111.96           100.0      25.0       2.1       1.0
"""
function get_forcing(simulation::DiscretizedSPAC) # returns forcing [..., ..., ..., ..., ...]

    # Extract data from simulation object
    # input forcings (from definition of parametrized SPAC):
    input_meteo = simulation.parametrizedSPAC.forcing.meteo
    # simulation.parametrizedSPAC.forcing.meteoiso
    simulation.parametrizedSPAC.forcing.storm_durations

    # processed forcings (integrated as time-varying parameters in ODEProblem):
    p_fT = simulation.ODESolution.prob.p;

    # handle time
    t_ref = p_fT.REFERENCE_DATE
    timepoints = simulation.parametrizedSPAC.forcing.meteo["p_days"]         # forcing period can be longer than simulation output (e.g. when using a spinup period)
    # timepoints = range(simulation.ODEProblem.tspan..., step = 1)           # Plot forcing as daily, even if solution output (ODESolution.t) is not dense
    # timepoints = range(extrema(simulation.ODESolution.t)..., step = 1)     # Plot forcing as daily, even if solution output (ODESolution.t) is not dense

    # 1) prepare data to plot
    # i) meteo forcing (wind, vappres, globrad, tmax, tmin, prec, lai, )
    return DataFrame(
        :dates                   => RelativeDaysFloat2DateTime.(timepoints, t_ref),
        :globrad_MJDayM2         => p_fT.p_GLOBRAD.(timepoints),
        :tmax_degC               => p_fT.p_TMAX.(timepoints),
        :tmin_degC               => p_fT.p_TMIN.(timepoints),
        :vappres_kPa             => p_fT.p_VAPPRES.(timepoints),
        :windspeed_ms            => p_fT.p_WIND.(timepoints),
        :prec_mmDay              => p_fT.p_PREC.(timepoints),
        :precdelta18O_permil     => p_fT.p_δ18O_PREC.(timepoints),
        :precdelta2H_permil      => p_fT.p_δ2H_PREC.(timepoints),
        :irrig_mmDay             => p_fT.p_IRRIG.(timepoints),
        :irrigdelta18O_permil    => p_fT.p_δ18O_IRRIG.(timepoints),
        :irrigdelta2H_permil     => p_fT.p_δ2H_IRRIG.(timepoints),
        :densef_percent          => p_fT.p_DENSEF.(timepoints)*100,
        :height_m                => p_fT.p_HEIGHT.(timepoints),
        :lai_m2m2                => p_fT.p_LAI.(timepoints),
        :sai_m2m2                => p_fT.p_SAI.(timepoints))

        # # ii) states (scalar and vector)
        # (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) = LWFBrook90.intern___get_SWATI_derivatives(simulation)
        # timepoints_states = simulation.ODESolution.t
        # # x3 = simulation.ODESolution_datetime;
        # x3 = RelativeDaysFloat2DateTime.(timepoints_states, t_ref);
        # y31 = hcat(reduce(hcat,  [simulation.ODESolution(t).INTR.mm   for t in simulation.ODESolution.t])',
        #             reduce(hcat, [simulation.ODESolution(t).INTS.mm   for t in simulation.ODESolution.t])',
        #             reduce(hcat, [simulation.ODESolution(t).SNOW.mm   for t in simulation.ODESolution.t])',
        #             reduce(hcat, [simulation.ODESolution(t).SNOWLQ.mm   for t in simulation.ODESolution.t])',
        #             reduce(hcat, [simulation.ODESolution(t).XYLEM.mm   for t in simulation.ODESolution.t])');
        # lbl31 = ["INTR [mm]" "INTS [mm]" "SNOW [mm]" "SNOWLQ [mm]" "XYLEM [mm]"];
        # y32 = hcat(sum(u_SWATI; dims=1)',
        #             reduce(hcat, [simulation.ODESolution(t).GWAT.mm   for t in simulation.ODESolution.t])');
        # lbl32 = ["Total Soil Water [mm]" "GWAT [mm]"];

        # # iii) fluxes
        # x4 = simulation.ODESolution_datetime;
        # y41 = hcat(  reduce(hcat, [simulation.ODESolution(t).RWU.mmday   for t in simulation.ODESolution.t])',
        #                                     sum(reduce(hcat, [simulation.ODESolution(t).TRANI.mmday   for t in simulation.ODESolution.t]); dims=1)')
        # lbl41 = ["RWU (mm/day)" "sum(TRANI)"]
end

"""
    get_soil_(symbols, simulation; depths_to_read_out_mm = nothing, days_to_read_out_d = nothing))

Returns a 2D DataFrame of soil variables with soil layers as columns and time steps as rows.
Supports a number of variables:
- `:θ`, `:theta` = volumetric soil moisture values (m3/m3)
- `:ψ`, `:psi` = soil matric potential (kPa)
- `:W` = soil wetness (-)
- `:SWATI` = soil water volumes contained in discretized layers (mm)
- `:K` = soil hydraulic conductivities (mm/day)
- `:δ18O`, `:δ2H`, `:d18O`, `:d2H` = isotopic signatures (delta)
- `:TRANI`, `:RWU` = root water uptake flux from each cell (mm/day)

The user can define timesteps as `days_to_read_out_d` or specific depths as `depths_to_read_out_mm`,
that are both optionally provided as numeric vectors, e.g. `depths_to_read_out_mm = [100, 150]` or `days_to_read_out_d = 1:1.0:100`
Function to read out soil variables from a simulated simulation.
- `symbols`: can be a single symbol (:θ) or a vector of symbols [:θ] or [:θ, :ψ,], see above for accepted symbols, and not that non-Unicode symbols are supported: e.g. [:theta, :psi, :delta18O, :delta2H, :d18O, :d2H]
- `simulation`: is a `DiscretizedSPAC` that has been simulated.
- `depths_to_read_out_mm`: either `nothing` or vector of Integers.
- `days_to_read_out_d`: either  `nothing` or vector of Floats representing days.

Examples

    get_soil_(:θ, simulation)
    get_soil_([:θ, :ψ, :K], simulation; depths_to_read_out_mm = [100, 200, 500, 1200])
"""
function get_soil_(
    symbols, simulation::DiscretizedSPAC;
    depths_to_read_out_mm = nothing,
    days_to_read_out_d = nothing,
    flag_return_Matrix = false) # legacy flag that can be set to request Matrix output... will be deprecated in future...

    solution          = simulation.ODESolution
    simulate_isotopes = simulation.parametrizedSPAC.solver_options.simulate_isotopes

    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"
    @assert (isnothing(depths_to_read_out_mm) || eltype(depths_to_read_out_mm) <: Integer) "`depths_to_read_out_mm` must be Vector{Int}. They are $(typeof(depths_to_read_out_mm))"
    @assert (isnothing(days_to_read_out_d) || all(days_to_read_out_d .>= 0)) "`days_to_read_out_d` must be Vector{Float64} and all >= 0. Received $(days_to_read_out_d)"

    # parse requested time points
    timepoints = isnothing(days_to_read_out_d) ? solution.t : days_to_read_out_d;
    timepoints = (timepoints isa AbstractArray ? timepoints : [timepoints]) # transform to vector even if input as scalar
    symbols    = (symbols    isa AbstractArray ? symbols    : [symbols]   ) # transform to vector even if input as scalar

    # 1a) get states
    # get auxiliary variables with requested time resolution (i.e. days_to_read_out_d)
    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        LWFBrook90.intern___get_SWATI_derivatives(simulation; days_to_read_out_d = days_to_read_out_d)
    # 1b) get isotope values
    if simulate_isotopes
        u_SWATI_d18O = reduce(hcat, [solution(t_days).SWATI.d18O for t_days = timepoints])
        u_SWATI_d2H  = reduce(hcat, [solution(t_days).SWATI.d2H  for t_days = timepoints])
    else
        u_SWATI_d18O = fill(missing, size(u_SWATI))
        u_SWATI_d2H  = fill(missing, size(u_SWATI))
        # old variant: # remove requested symbols from
        # old variant: symbols[symbols .!= (:δ18O) .&& symbols .!= (:delta18O) .&& symbols .!= (:delta2H) .&& symbols .!= (:delta2H)]
    end
    # 2) get fluxes
    u_TRANI = reduce(hcat, [solution(t_days).TRANI.mmday for t_days = timepoints])
    u_PLRFI = reduce(hcat, [solution(t_days).PLRFI.mmday for t_days = timepoints])
    # TODO: possibly add BYFLI, INFLI, DSFLI

    # Setup DataFrame to fill
    df = DataFrame()
    df[:, :time] = timepoints

    # Fill DataFrame with requested soil layers and requested variables
    if isnothing(depths_to_read_out_mm)
        lower_boundaries = round.(Int, cumsum(simulation.ODESolution.prob.p.p_soil.p_THICK))      # NOTE: this rounds hardcoded to mm
        # center = round.(Int, lower_boundaries - 0.5*simulation.ODESolution.prob.p.p_soil.p_THICK) # NOTE: this rounds hardcoded to mm
        depths_to_read_out_mm = lower_boundaries
    end
    idxs = LWFBrook90.intern___get_soil_idx(simulation, depths_to_read_out_mm; only_valid_idxs = true)

    for symbol in sort(symbols)
        variable_to_return, unit_to_return = (
            symbol in [:θ, :theta              ] ? (u_aux_θ,      "m3m3"   ) :   # also accept non-unicode
            symbol in [:ψ, :psi                ] ? (u_aux_PSIM,   "kPa"    ) :   # also accept non-unicode
            # symbol == :ψtot ? u_aux_PSITI :
            symbol in [:δ18O, :delta18O, :d18O ] ? (u_SWATI_d18O, "permil" ) : # also accept non-unicode
            symbol in [:δ2H,  :delta2H,  :d2H  ] ? (u_SWATI_d2H,  "permil" ) : # also accept non-unicode
            symbol in [:W                      ] ? (u_aux_WETNES, "_"      ) :
            symbol in [:SWATI                  ] ? (u_SWATI,      "mm"     ) :
            symbol in [:K                      ] ? (p_fu_KK,      "mmday"  ) :
            symbol in [:TRANI, :RWU            ] ? (u_TRANI,      "mmday"  ) :
            symbol in [:PLRFI                  ] ? (u_PLRFI,      "mmday"  ) :
            # TODO: possibly add BYFLI, INFLI, DSFLI
            error("Unknown symobl $(symbol) requested in `get_soil_()`"))

        # Add columns
        # append requested layers with the depths as column titles
        [df[:, Symbol("$(string(symbol))_$(unit_to_return)_$(Int(k))mm")] = variable_to_return[v,:] for
            (k,v) in sort(collect(idxs)) if v != 0]; # Note we sort the idxs by depth
    end

    # Return as DataFrame or Matrix (latter is to support legacy code)
    if (flag_return_Matrix)
        @assert length(symbols)==1 """
            Aborting `get_soil_()` with non-DataFrame output for more than 1 soil variable. It is ambiguous which colum represents which variable (θ, ψ, ...).
            Correct by either requesting only 1 variable type by supplying only one `symbols`. Alternatively mutliple variable types are supported for specified `depths_to_read_out_mm`.
        """
        return permutedims(Matrix(df[:,Not(:time)]))
    else
        return df
    end
end

"""
    get_water_partitioning(simulation)

Returns three 2D DataFrame of water fluxes with different fluxes as columns and time steps as rows.
The three DataFrames are in daily, monthly and yearly resolution and span the entire simulation.

Examples
    df_part_daily, df_part_monthly, df_part_yearly = get_water_partitioning(simulation)
"""
function get_water_partitioning(simulation::DiscretizedSPAC;)
    solution          = simulation.ODESolution
    simulate_isotopes = simulation.parametrizedSPAC.solver_options.simulate_isotopes

    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"
    @assert all(diff(simulation.ODESolution_datetime) .== Millisecond(Day(1))) """
    Solution is not computed with daily output resolution.
    Make sure you provide simulate!(save_everystep = false, saveat = ...) with `saveat` in daily resolution."""

    ks = keys(simulation.ODESolution.u[1].accum)
    df_partitioning_raw =
        DataFrame((:date => simulation.ODESolution_datetime ),
                (k => [ut.accum[k] for ut in simulation.ODESolution.u] for k in ks)...)
    # Compute ETa, Es, Esn, Ei, Ta, P, Td, D, R, Swat
    df_partitioning_daily = @chain df_partitioning_raw begin
        @rtransform begin
            :ETa           = :evap # is actual evapotranspiration, i.e. sum of IRVP + ISVP + SNVP + SLVP + sum(aux_du_TRANI)
            :Esoil         = :cum_d_slvp
            :Esnow         = :cum_d_snvp
            :Einterception = :cum_d_irvp + :cum_d_isvp
            :Ta            = :cum_d_tran
            :Precip        = :cum_d_prec
            :Irrig         = :cum_d_irrig
            :Td            = :cum_d_ptran - :cum_d_tran
            :D             = - :vrfln
            # :R1          = -(:flow - :vrfln)    # flow = srfl+byfl+dsfli+gwfl, gwfl, vrfln
            :R             = -(:srfl + :byfl + :dsfl) # This is more correctly not accounting for gwfl, thereby excluding state variable GWAT from the balance
            :Swat          = :StorageSWAT
        end
        @transform :year = year.(:date)
        @transform :month = month.(:date)
        @select(:date, :year, :month, :ETa,:Esoil,:Esnow,:Einterception,:Ta,:Precip,:Irrig,:Td,:D,:R,:Swat)
    end

    # Aggregate to monghly and yearly
    df_partitioning_monthly = @chain df_partitioning_daily begin
        groupby([:year, :month])
        combine(
            nrow,
            [:ETa,:Esoil,:Esnow,:Einterception,:Ta,:Precip,:Irrig,:Td,:D,:R] .=> sum,
            [:Swat] .=> mean, renamecols=false)
        @rtransform :date = Date(:year, :month)
        select(Between(:year, :nrow), :date, All()) # Bring date to beginning
    end
    df_partitioning_yearly = @chain df_partitioning_daily begin
        groupby([:year])
        combine(
            nrow,
            [:ETa,:Esoil,:Esnow,:Einterception,:Ta,:Precip,:Irrig,:Td,:D,:R] .=> sum,
            [:Swat] .=> mean, renamecols=false)
        @rtransform :date = Date(:year)
        select(Between(:year, :nrow), :date, All()) # Bring date to beginning
    end
    # and also define color palette
    return (df_partitioning_daily, df_partitioning_monthly, df_partitioning_yearly)
end
#######
####### Internal API =======================================================================
#######

function intern___get_soil_idx(simulation::DiscretizedSPAC, depths_to_read_out_mm; only_valid_idxs = false)
    @assert all(depths_to_read_out_mm .> 0) # depths and lower_boundaries must all be positive numbers
    lower_boundaries = round.(Int, cumsum(simulation.ODESolution.prob.p.p_soil.p_THICK))      # NOTE: this rounds hardcoded to mm
    # center = round.(Int, lower_boundaries - 0.5*simulation.ODESolution.prob.p.p_soil.p_THICK) # NOTE: this rounds hardcoded to mm

    idx_to_read_out = fill(0, length(depths_to_read_out_mm))
    for (it, curr_depth_mm) in enumerate(depths_to_read_out_mm)
        if (curr_depth_mm > maximum(lower_boundaries))
            # Only read out values that are within the simulation domain
            idx_to_read_out[it] = 0 # 0 means this depth has not been simulated
            @warn "Requested read-out depth of $curr_depth_mm is below simulation domain and is silently omitted for the output."
        else
            idx_to_read_out[it] = findfirst(curr_depth_mm .<= lower_boundaries)
        end
    end
    all_idxs = Dict((d => i) for (d,i) in zip(depths_to_read_out_mm, idx_to_read_out))
    if only_valid_idxs
        return filter((k, v)::Pair -> v != 0, all_idxs)
    else # default
        return all_idxs
    end
end
function intern___get_SWATI_derivatives(simulation::DiscretizedSPAC; days_to_read_out_d = nothing)
    solution = simulation.ODESolution
    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"

    p_soil = solution.prob.p.p_soil
    NLAYER = p_soil.NLAYER
    if isnothing(days_to_read_out_d)
        u_SWATI = reduce(hcat, [solution.u[t_idx].SWATI.mm  for t_idx = eachindex(solution)])
    else
        u_SWATI = reduce(hcat, [solution(t_days).SWATI.mm for t_days = days_to_read_out_d])
    end

    u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK =
        (fill(NaN, size(u_SWATI)) for i in 1:5)
    for t in 1:size(u_SWATI,2)
        (u_aux_WETNES[:, t], u_aux_PSIM[:, t], u_aux_PSITI[:, t], u_aux_θ[:, t], p_fu_KK[:, t]) =
            KPT.derive_auxiliary_SOILVAR(u_SWATI[:,t], p_soil)
    end
    # returns arrays of dimenstion (t,z) where t is number of timesteps and z number of computational layers
    return (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK)
end
function intern___get_scalars(compartments_to_extract, units_to_extract, simulation::DiscretizedSPAC, days_to_read_out_d = nothing)
    solution = simulation.ODESolution;
    @assert !isnothing(solution) "Solution was not yet computed. Please simulate!(simulation)"
    @assert (isnothing(days_to_read_out_d) || all(days_to_read_out_d .>= 0)) "`days_to_read_out_d` must be Vector{Float64} and all >= 0. Received $(days_to_read_out_d)"

    # parse requested time points and compartments
    timepoints = isnothing(days_to_read_out_d) ? solution.t : days_to_read_out_d;
    df_time = DataFrame(:time => (timepoints isa AbstractArray ? timepoints : [timepoints])) # transform to vector even if input as scalar

    # Extract scalar values from state variable vector u
    # Treat precipitation separately
    is_prec = compartments_to_extract .== :PREC
    cycle_over1 = zip(compartments_to_extract[Not(is_prec)], units_to_extract[Not(is_prec)])
    cycle_over2 = zip(compartments_to_extract[is_prec], units_to_extract[is_prec])

    # u_aboveground = collect(solution(t_days)[comp][uni]  for t_days = days_to_read_out_d, (comp,uni) in cycle_over1)
    df_states = DataFrame((
        string(comp) * "_" * string(uni) => [
            solution(t_days)[comp][uni] for t_days = timepoints
        ] for (comp, uni) in cycle_over1)...)

    # Extract precipitation from forcing
    df_PREC = DataFrame((string(comp) * "_" * string(uni) =>
        # comp == :PREC ?
        uni ∈ [:mmday,  :mm] ? solution.prob.p.p_PREC.(timepoints) :
            uni ∈ [:d18O, :δ18O] ? solution.prob.p.p_δ18O_PREC.(timepoints) :
            uni ∈ [:d2H , :δ2H ] ? solution.prob.p.p_δ2H_PREC.(timepoints) :
            # uni ∈ [:mmday,  :mm] ? solution.prob.p.p_IRRIG.(timepoints) :
            # uni ∈ [:d18O, :δ18O] ? solution.prob.p.p_δ18O_IRRIG.(timepoints) :
            # uni ∈ [:d2H , :δ2H ] ? solution.prob.p.p_δ2H_IRRIG.(timepoints) :
            error("Unknown compartments or units to extract provided.")
        for (comp, uni) in cycle_over2)...)

    # TODO: if we cycle only once over states and fluxes we would get the same order as provided in
    #       compartments_to_extract, units_to_extract. This would be better behavior.
    #       But only modify this if intern___get_scalars() becomes part of the public API

    # returns DataFrame of dimenstion (t, N_variables = 9)
    #   where t is number of time steps
    #   and 9 = 1 (time) + 1(PREC) + 6 (number of compartments)
    return hcat(df_time, df_PREC, df_states)
end
function intern___get_water_partitioning_colorpalette()
        color_palette_Meusburger2022 = reverse([
            "Td"            => :red2,
            "Ta"            => :darkolivegreen2,
            "Einterception" => :forestgreen,
            "Esoil"         => :khaki3,
            "Esnow"         => :white,
            "R"             => :lightskyblue,
            "D"             => :steelblue4,
            # "ETa"   => :black,
            # "P2"    => :darkblue,
            "P"     => :darkblue,
            # "Swat" => :brown
        ])
        color_palette_Schmidt_Walter2020 = reverse([
            "Td"            => colorant"#d01c8b", #:palevioletred4,
            "Ta"            => colorant"#abdda4", #:darkseagreen,
            "Einterception" => colorant"#fdae61", #:lightyellow,
            "Esoil"         => colorant"#FFD700", #:navajowhite, #:orange2,
            "Esnow"         => :white,
            "R"             => colorant"#91bfdb", #:slategray2,
            "D"             => colorant"#2b83ba", #:skyblue4,
            "Precip"        => :black,
            # "ETa"  => :black,
            # "P2"   => :darkblue,
            # "Swat" => :brown
        ])
        return color_palette = color_palette_Schmidt_Walter2020
end
function intern___get_RWU_centroid(rows_RWU_mmDay, y_center)
    # old solution, that gives bad values when some RWU is negative:
    RWU_percent = rows_RWU_mmDay ./ sum(rows_RWU_mmDay; dims = 1)
    RWUcentroidLabel = "mean RWU depth"
    if (any(RWU_percent .< 0))
        @warn "Some root water outfluxes detected. Centroid of RWU is  based only on uptakes."
        rows_RWU_mmDay_onlyUptake = ifelse.(rows_RWU_mmDay.>0,rows_RWU_mmDay, 0)
        RWU_percent_onlyUptake = rows_RWU_mmDay_onlyUptake ./ sum(rows_RWU_mmDay_onlyUptake; dims = 1)
        RWU_percent = RWU_percent_onlyUptake
        # RWU_percent = min.(0,rows_RWU_mmDay) ./ sum(min.(0,rows_RWU_mmDay); dims = 1)

        RWUcentroidLabel = "mean RWU depth\n(based on uptake only)"
    end

    return (row_RWU_centroid_mm = sum(RWU_percent .* y_center; dims=1),
            RWUcentroidLabel)
end

function intern___get_data_for_isotopePlot(simulation)
    # 1) prepare data to plot
    solu = simulation.ODESolution

    # 1a) extract data from solution object `solu`
    # # Two ways to extract data from soil object: using `[]` or `()`
    # u_SWATI = reduce(hcat, [solu[t_idx].SWATI.mm  for t_idx = eachindex(solu)])
    # u_SWATI = reduce(hcat, [solu(t_days).SWATI.mm for t_days = days_to_read_out_d])

    # days_to_read_out_d decides which points to use:
    # days_to_read_out_d = nothing # read out all simulation steps
    days_to_read_out_d = :daily  # read out only daily values
    if isnothing(days_to_read_out_d)
        days_to_read_out_d = solu.t # warning this can
        error("Too many output times requested. This easily freezes the program...")
    elseif days_to_read_out_d == :daily
        days_to_read_out_d = unique(round.(solu.t))
    end
    t_ref = solu.prob.p.REFERENCE_DATE
    x = RelativeDaysFloat2DateTime.(days_to_read_out_d, t_ref);
    y_center = cumsum(solu.prob.p.p_soil.p_THICK) - solu.prob.p.p_soil.p_THICK/2

    simulate_isotopes = solu.prob.p.simulate_isotopes
    @assert simulate_isotopes "Provided DiscretizedSPAC() did not simulate isotopes"

    # # Some hardcoded options:
    # xlimits = RelativeDaysFloat2DateTime.(solu.prob.tspan, t_ref)
    # tick_function = (x1, x2) -> PlotUtils.optimize_ticks(x1, x2; k_min = 4)
    # # color_scheme = :default # https://docs.juliaplots.org/latest/generated/colorschemes
    # # color_scheme = :blues
    # color_scheme = :heat

    # Results
    t_ref = solu.prob.p.REFERENCE_DATE
    t1 = range(extrema(solu.t)..., step = 1) # Plot forcing as daily, even if solution output (ODESolution.t) is not dense
    # x1 = RelativeDaysFloat2DateTime.(t1, t_ref);
    # maxdepth = maximum(cumsum(solu.prob.p.p_soil.p_THICK))
    # Scalar quantities (use DataFrame for Storage, reshape to Matrix/Rows for plotting)
    df_scalar = DataFrame(
        days_to_read_out_d = days_to_read_out_d,
        col_PREC_amt_dense  = solu.prob.p.p_PREC.(t1),
        col_PREC_d18O_dense = solu.prob.p.p_δ18O_PREC.(t1),
        col_PREC_d2H_dense  = solu.prob.p.p_δ2H_PREC.(t1),
        col_PREC_d18O = solu.prob.p.p_δ18O_PREC.(days_to_read_out_d),
        col_PREC_d2H  = solu.prob.p.p_δ2H_PREC.(days_to_read_out_d),
        col_IRRIG_amt_dense  = solu.prob.p.p_IRRIG.(t1),
        col_IRRIG_d18O_dense = solu.prob.p.p_δ18O_IRRIG.(t1),           # TODO: do we really want to output nothing? Or simply NaN? (Prevents Union type: Union{Nothing, FLoat64})
        col_IRRIG_d2H_dense  = solu.prob.p.p_δ2H_IRRIG.(t1),            # TODO: do we really want to output nothing? Or simply NaN? (Prevents Union type: Union{Nothing, FLoat64})
        col_IRRIG_d18O = solu.prob.p.p_δ18O_IRRIG.(days_to_read_out_d), # TODO: do we really want to output nothing? Or simply NaN? (Prevents Union type: Union{Nothing, FLoat64})
        col_IRRIG_d2H  = solu.prob.p.p_δ2H_IRRIG.(days_to_read_out_d),  # TODO: do we really want to output nothing? Or simply NaN? (Prevents Union type: Union{Nothing, FLoat64})
        col_INTS_d18O = [solu(t).INTS.d18O for t in days_to_read_out_d],
        col_INTR_d18O = [solu(t).INTR.d18O for t in days_to_read_out_d],
        col_SNOW_d18O = [solu(t).SNOW.d18O for t in days_to_read_out_d],
        col_GWAT_d18O = [solu(t).GWAT.d18O for t in days_to_read_out_d],
        col_RWU_d18O  = [solu(t).RWU.d18O for t in days_to_read_out_d],
        col_XYL_d18O  = [solu(t).XYLEM.d18O for t in days_to_read_out_d],
        col_INTS_d2H  = [solu(t).INTS.d2H for t in days_to_read_out_d],
        col_INTR_d2H  = [solu(t).INTR.d2H for t in days_to_read_out_d],
        col_SNOW_d2H  = [solu(t).SNOW.d2H for t in days_to_read_out_d],
        col_GWAT_d2H  = [solu(t).GWAT.d2H for t in days_to_read_out_d],
        col_RWU_d2H   = [solu(t).RWU.d2H for t in days_to_read_out_d],
        col_XYL_d2H   = [solu(t).XYLEM.d2H for t in days_to_read_out_d])

    # Vector quantities
    cols_SWAT_d18O = reduce(vcat, [solu(t).SWATI.d18O' for t in days_to_read_out_d])
    cols_SWAT_d2H  = reduce(vcat, [solu(t).SWATI.d2H' for t in days_to_read_out_d])

    # Compute RWU centroid
    cols_RWU_mmDay  = reduce(vcat, [solu(t).TRANI.mmday'   for t in days_to_read_out_d])
    rows_RWU_mmDay  = permutedims(cols_RWU_mmDay)
    row_RWU_centroid_mm, RWUcentroidLabel = LWFBrook90.intern___get_RWU_centroid(rows_RWU_mmDay, y_center)
    col_RWU_centroid_mm = reshape(row_RWU_centroid_mm, :)
    # if (RWUcentroid == :showRWUcentroid)
    # end

    # Finalize DataFrame for storage
    df_scalar.col_RWU_centroid_mm = col_RWU_centroid_mm
    # df_amt  = DataFrame(cols_SWAT_amt,  ["amt_$depth"  for depth in (round.(Int, y_center))])
    df_d18O = DataFrame(cols_SWAT_d18O, ["d18O_$depth" for depth in (round.(Int, y_center))])
    df_d2H  = DataFrame(cols_SWAT_d2H,  ["d2H_$depth"  for depth in (round.(Int, y_center))])
    df = hcat(df_scalar, df_d18O, df_d2H)

    return (df = df, RWUcentroidLabel = RWUcentroidLabel)
end

#######
####### Public facing API: plotting recipes ================================================
#######
function plot_monthly_water_partitioning(df_partitioning_monthly, fig = Figure())
    color_palette = intern___get_water_partitioning_colorpalette()
    # Preprocess
    df_partitioning_monthly_forMakie = @chain df_partitioning_monthly begin
        stack(Not([:year, :month, :nrow, :date]))
        # only keep variables we need
        @subset(:variable .∈ ([first(pair) for pair in color_palette],))
        # make categorical
        @transform :variable = CategoricalArrays.categorical(:variable, levels = [first(pair) for pair in color_palette])
        @transform :variable_code = CategoricalArrays.levelcode.(:variable)
        # Remove fluxes that were not computed (e.g. removes runoff)
        @subset(:value .!= 0.0)
        end
    # Plot
    aog_monthly = AlgebraOfGraphics.mapping(
            :date => "",
            :value => "Water flux per month (mm)",
            stack = :variable,
            color = :variable => "") *
        (# bar plot of fluxes
        AlgebraOfGraphics.data(@subset(df_partitioning_monthly_forMakie, :variable .!= "Precip")) * AlgebraOfGraphics.visual(BarPlot, gap = -31*0.8) +
        # line plot of precip input
        AlgebraOfGraphics.data(@subset(df_partitioning_monthly_forMakie, :variable .== "Precip")) * AlgebraOfGraphics.visual(Lines)
        )
    xticks = sort(unique(Dates.floor.(df_partitioning_monthly_forMakie.date, Dates.Month(6))))

    aog_draw = AlgebraOfGraphics.draw!(fig, aog_monthly, palettes = (; color = color_palette),
        # axis = (; xticks = AlgebraOfGraphics.datetimeticks((x -> Dates.format(x, "mm\nY")), (Date.(xticks)))))
        axis = (; ygridvisible = true,
                xticks = AlgebraOfGraphics.datetimeticks((x -> Dates.format(x, "u\nY")), (Date.(xticks)))))
    return fig, aog_draw
end
function plot_yearly_water_partitioning(df_partitioning_yearly, fig = Figure())
    color_palette = intern___get_water_partitioning_colorpalette()
    # Preprocess
    df_partitioning_yearly_forMakie = @chain df_partitioning_yearly begin
        stack(Not([:year, :nrow, :date]))
        # only keep variables we need
        @subset(:variable .∈ ([first(pair) for pair in color_palette],))
        # make categorical
        @transform :variable = CategoricalArrays.categorical(:variable, levels = [first(pair) for pair in color_palette])
        @transform :variable_code = CategoricalArrays.levelcode.(:variable)
        # Remove fluxes that were not computed (e.g. removes runoff)
        @subset(:value .!= 0.0)
        end
    # Plot
    aog_yearly = AlgebraOfGraphics.mapping(
        :date => "",
        :value => "Water flux per year (mm)",
        stack = :variable,
        color = :variable => "") *
    (# bar plot of fluxes
    AlgebraOfGraphics.data(@subset(df_partitioning_yearly_forMakie, :variable .!= "Precip")) * AlgebraOfGraphics.visual(BarPlot, gap = -366*0.8) +
    # line plot of precip input
    AlgebraOfGraphics.data(@subset(df_partitioning_yearly_forMakie, :variable .== "Precip")) * AlgebraOfGraphics.visual(Lines)
    )

    xticks = sort(unique(df_partitioning_yearly_forMakie.year))
    aog_draw = AlgebraOfGraphics.draw!(fig, aog_yearly,
        palettes = (; color = color_palette),
        axis = (; ygridvisible = true,
                #   xticks = AlgebraOfGraphics.datetimeticks((x -> Dates.format(x, "Y-mm")), (Date.(xticks)))))
                  xticks = AlgebraOfGraphics.datetimeticks((x -> Dates.format(x, "Y")), (Date.(xticks)))))
    return fig, aog_draw
end















"""
    plotamounts(simulation::DiscretizedSPAC)
    plotamounts(simulation::DiscretizedSPAC, compartments::Symbol)
    plotamounts(simulation::DiscretizedSPAC, compartments::Symbol, RWUcentroid::Symbol)

Plots the amount results of a SPAC Simulation. By default both above and belowground.
The user can override this with the second argument isotope as one of `:aboveground`, `:belowground`, or `:above_and_belowground`.
RWUcentroid can have values of either `:dontShowRWUcentroid` or `:showRWUcentroid`.
"""
@userplot PlotAmounts
@recipe function f(plam::PlotAmounts)
    # 0) parse input arguments

    if length(plam.args) == 3
        simulation = plam.args[1]
        compartments = plam.args[2]
        RWUcentroid = plam.args[3]
    elseif length(plam.args) == 2
        simulation = plam.args[1]
        compartments = plam.args[2]
        RWUcentroid = :dontShowRWUcentroid
    elseif length(plam.args) == 1
        simulation = plam.args[1]
        compartments = :above_and_belowground
        RWUcentroid= :dontShowRWUcentroid
    else
        error("plotamounts requires an unnamed first argument of type DiscretizedSPAC, and optional unnamed second/third arguments (:aboveground, :belowground, or :above_and_belowground) and (:dontShowRWUcentroid, :showRWUcentroid). Other arguments to plot() should be separated by `;`.")
    end
    if (compartments != :above_and_belowground)
        error("TODO: currently selecting a subset of compartments are not supported by plotamounts. Please only do: `plotamounts(simulation)` or `plotamounts(simulation, :above_and_belowground, :showRWUcentroid)`")
    end
    if !(RWUcentroid == :dontShowRWUcentroid || RWUcentroid == :showRWUcentroid)
        error("Third unnamed argument to plotamounts should be one of (:dontShowRWUcentroid, :showRWUcentroid). Got: $(RWUcentroid)")
    end
    if !(compartments == :aboveground || compartments == :belowground || compartments == :above_and_belowground)
        error("Second unnamed argument to plotamounts should be one of (:above_and_belowground, :aboveground, or :belowground). Got: $( compartments)")
    end
    if !(simulation isa DiscretizedSPAC)
        error("First unnamed argument to plotamounts should be of type DiscretizedSPAC. Got: $(typeof(simulation))")
    end
    if isnothing(simulation.ODESolution)
        error("plotamounts requires a solved system. Please `simulate!()` the DiscretizedSPAC first.")
    end

    # 1) prepare data to plot
    solu = simulation.ODESolution

    # 1a) extract data from solution object `solu`
        # # Two ways to extract data from soil object: using `[]` or `()`
        # u_SWATI = reduce(hcat, [solu[t_idx].SWATI.mm  for t_idx = eachindex(solu)])
        # u_SWATI = reduce(hcat, [solu(t_days).SWATI.mm for t_days = days_to_read_out_d])

        # days_to_read_out_d decides which points to use:
        # days_to_read_out_d = nothing # read out all simulation steps
        days_to_read_out_d = :daily  # read out only daily values
        if isnothing(days_to_read_out_d)
            days_to_read_out_d = solu.t # warning this can
            error("Too many output times requested. This easily freezes the program...")
        elseif days_to_read_out_d == :daily
            days_to_read_out_d = unique(round.(solu.t))
        end
        t_ref = solu.prob.p.REFERENCE_DATE
        x = RelativeDaysFloat2DateTime.(days_to_read_out_d, t_ref);
        y_center = cumsum(solu.prob.p.p_soil.p_THICK) - solu.prob.p.p_soil.p_THICK/2

    # Some hardcoded options:
    xlimits = RelativeDaysFloat2DateTime.(solu.prob.tspan, t_ref)
    # color_scheme = :default # https://docs.juliaplots.org/latest/generated/colorschemes
    # color_scheme = :heat
    color_scheme = :blues
    """
        cgrad_two_separate_gradients_at_mid()
    """
    function cgrad_two_separate_gradients_at_mid(z, mid;
        # colors = [:white, :blue, :red, :black],
        colors = [colorant"lightblue", colorant"darkblue", :black, :red],
        up     = maximum(z),
        down   = minimum(z),
        kwargs...)

        eps_ = eps()*abs(up - down)        # scale difference that it remains tractable
        # if mid >= up; mid = up - eps_; end # if red part not in data, correct accordingly
        if mid >= up;
            # if red part not in data, use gradient only with blue
            return cgrad(colors[1:2], [0, 1]; kwargs...)
        else
            @assert down <= mid - eps_ < mid <= up
            values = ([down, mid - eps_, mid, up] .- down) ./ (up - down) # scale to range 0-1
            # @show colors
            # @show values
            cgrad(colors, values; kwargs...)
        end
    end
            # Not needed if we're moving to Makie...
            # function log_modulus(x; base = exp(1))
            #     sign.(x) .* log.(base,  abs.(x) .+ 1)
            # end
            # function invlog_modulus(x; base = exp(1))
            #     sign.(x) .* (base .^ (abs.(x)) .- 1)
            # end
            # # rows_ψₘpF_logmodulus  = log10.(u_aux_PSIM) #(from kPa to log10(hPa))
            # #     # from: https://discourse.julialang.org/t/is-there-a-way-to-plot-negative-numbers-on-a-log-scale/49674/7?u=fabern
            # #     # as long as it is not yet implemented: https://github.com/JuliaPlots/Plots.jl/issues/1235
            # # z = sign.(y).*log10.(abs.(y).+1)
            # # plot(x, z, color = z,
            # #     yformatter =  x -> "$(sign(x) < 0 ? "-" : "")10^$x")

    true_to_check_colorbar = true; # set this flag to false for final plot, true for debugging.
    tick_function = (x1, x2) -> PlotUtils.optimize_ticks(x1, x2; k_min = 4)
    y_soil_ticks = tick_function(0., round(maximum(cumsum(solu.prob.p.p_soil.p_THICK))))[1]
    y_ticks    = [-500;       -300;       -200;       -100;          y_soil_ticks;             (maximum(cumsum(solu.prob.p.p_soil.p_THICK)) .+ [    100;      250;     400])]
    y_labels   = ["PREC";   "INTS";     "INTR";     "SNOW";     round.(y_soil_ticks; digits=0);                                                "GWAT";    "RWU";     "XYLEM"]


    # Results
    t_ref = solu.prob.p.REFERENCE_DATE
    t1 = range(extrema(solu.t)..., step = 1) # Plot forcing as daily, even if solution output (ODESolution.t) is not dense
    x1 = RelativeDaysFloat2DateTime.(t1, t_ref);

    #t1_forcing = simulation.parametrizedSPAC.forcing.meteo["p_days"]
    #x1_input = RelativeDaysFloat2DateTime.(t1_forcing, t_ref);
    #row_PREC_amt_dense_input = reshape(simulation.parametrizedSPAC.forcing.meteo["p_PREC"](t1_forcing),1,:)
    # bar(reshape(x1_input,1,:), reshape(row_PREC_amt_dense_input,1,:))

    row_PREC_amt_dense = reshape(solu.prob.p.p_PREC.(t1), 1, :)
    # row_IRRIG_amt_dense = reshape(solu.prob.p.p_IRRIG.(t1), 1, :)
    rows_SWAT_amt  = reduce(hcat, [solu(t).SWATI.mm    for t in days_to_read_out_d])
    rows_RWU_mmDay = reduce(hcat, [solu(t).TRANI.mmday for t in days_to_read_out_d])
    row_NaN       = fill(NaN, 1,length(x))
    col_INTS_amt = [solu(t).INTS.mm for t in days_to_read_out_d]
    col_INTR_amt = [solu(t).INTR.mm   for t in days_to_read_out_d]
    col_SNOW_amt = [solu(t).SNOW.mm   for t in days_to_read_out_d]
    col_GWAT_amt = [solu(t).GWAT.mm   for t in days_to_read_out_d]
    col_RWU_mmDay  = [solu(t).RWU.mmday for t in days_to_read_out_d]
    col_XYL_amt  = [solu(t).XYLEM.mm  for t in days_to_read_out_d]

    # For water balance errors
    # WB error defined by Ireson et al. 2023
    function compute_WB_Ireson2023(solution)
        # Compute In/Out Fluxes and Storages
        days = unique(round.(solution.t)) # we have to read this out in a regular grid!
        dat_daily_cumulativeFluxes_and_Storage = DataFrame(
            # NOTE: we need not to take differences as these cumulative fluxes were set to 0
            #       every day by a solver callback
            slfl  = [solu(t).accum.slfl for t in days],
            byfl  = [solu(t).accum.byfl for t in days],
            vrfln = [solu(t).accum.vrfln for t in days],
            dsfl  = [solu(t).accum.dsfl for t in days],
            swat  = vcat([sum(solu(t).SWATI.mm, dims=1) for t in days]...), # storage
            tran  = [solu(t).accum.cum_d_tran for t in days],
            slvp  = [solu(t).accum.cum_d_slvp for t in days],
            prec  = [solu(t).accum.cum_d_prec for t in days],
            irrig = [solu(t).accum.cum_d_irrig for t in days],
            evap  = [solu(t).accum.evap for t in days], # evap is sum of irvp, isvp, snvp, slvp, sum(trani)
            flow  = [solu(t).accum.flow for t in days], # flow is sum of byfli, dsfli, gwfl
            seep  = [solu(t).accum.seep for t in days],
            gwat  = [solu(t).GWAT.mm for t in days], # storage
            snow  = [solu(t).SNOW.mm for t in days], # storage
            intr  = [solu(t).INTR.mm for t in days], # storage
            ints  = [solu(t).INTS.mm for t in days], # storage
            )

        function compute_error_SWATI(df) #
            cumInflow_mm  = cumsum(df.slfl - df.byfl) # =INFL,               corresponds to q(t,0)  in Ireson 2023 eq 16 with additionally sources and sinks
            cumOutflow_mm = cumsum(df.vrfln + df.dsfl + df.tran + df.slvp) # corresponds to q(t,zN) in Ireson 2023 eq 16 with additionally sources and sinks
            error_mm = (cumInflow_mm .- cumOutflow_mm) .- (df.swat .- df.swat[1])
            return cumInflow_mm, cumOutflow_mm, df.swat, error_mm
        end
        function compute_error_ModelDomain(df) #
            cumInflow_mm  = cumsum(df.prec + df.irrig)          # corresponds to q(t,0)  in Ireson 2023 eq 16 with additionally sources and sinks
            cumOutflow_mm = cumsum(df.evap + df.flow + df.seep) # corresponds to q(t,zN) in Ireson 2023 eq 16 with additionally sources and sinks
            storage = df.swat .+ df.gwat + df.snow + df.intr + df.ints
            error_mm = (cumInflow_mm .- cumOutflow_mm) .- (storage .- storage[1])
            return cumInflow_mm, cumOutflow_mm, storage, error_mm
        end

        # Compute Water balance error
        # a) Cumulative error metric (time evolution of water balance error)
        cum_qIn_SWATI, cum_qOut_SWATI, storage_SWATI, εB_cumulative_SWATI = compute_error_SWATI(dat_daily_cumulativeFluxes_and_Storage)
        cum_qIn_ALL, cum_qOut_ALL, storage_ALL, εB_cumulative_ALL = compute_error_ModelDomain(dat_daily_cumulativeFluxes_and_Storage)

        # b) Scalar error metric
        εB_SWATI = εB_cumulative_SWATI[end]
        εB_ALL   = εB_cumulative_ALL[end]
        εR_SWATI = (diff(cum_qIn_SWATI) .- diff(cum_qOut_SWATI)) .- (diff(storage_SWATI)) # mm, Bias error for intervals t = (0,tM), Ireson 2023 eq. 15
        εR_SWATI = sqrt(sum(εR_SWATI.^2)/length(εR_SWATI))
        εR_ALL   = (diff(cum_qIn_ALL)   .- diff(cum_qOut_ALL))   .- (diff(storage_ALL))   # mm, Bias error for intervals t = (0,tM), Ireson 2023 eq. 15
        εR_ALL = sqrt(sum(εR_ALL.^2)/length(εR_ALL))

        return (εB_cumulative_SWATI, εB_cumulative_ALL, εB_SWATI, εB_ALL, εR_SWATI, εR_ALL,
                    (cum_qIn_SWATI, cum_qOut_SWATI, storage_SWATI, cum_qIn_ALL, cum_qOut_ALL, storage_ALL))
    end

    (εB_cumulative_SWATI, εB_cumulative_ALL, εB_SWATI, εB_ALL, εR_SWATI, εR_ALL,
        (cum_qIn_SWATI, cum_qOut_SWATI, SWATI_total, cum_qIn_ALL, cum_qOut_ALL, ALL_total)) =
        compute_WB_Ireson2023(solu);

    if (RWUcentroid == :showRWUcentroid)
        row_RWU_centroid_mm, RWUcentroidLabel = intern___get_RWU_centroid(rows_RWU_mmDay, y_center)
    end

    # reduce how deep to plot soil:
    soil_discr_to_plot = simulation.parametrizedSPAC.soil_discretization.df#[simulation.soil_discretization.Lower_m .>= -0.2, :]

    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        intern___get_SWATI_derivatives(simulation, days_to_read_out_d = days_to_read_out_d);

    # compute total amount of soil water: a) in whole domain, b) per physical soil layer
    # a)
    col_sumSWATI_amt = sum(u_SWATI; dims = [1])[:]

    # b)
    # sum up SWATI to get total SWAT for each of the physical soil horizons
    # terminology: soil horizon refers to a physical horizon, discretization refers to the computational layers/cells
    # soil_horizons = unique(soil_discr_to_plot, :HorizonNr)
    soil_horizons = unique(soil_discr_to_plot, :HorizonNr)[:,[:HorizonNr, :Horizon_Upper_m, :Horizon_Lower_m]]
    cols_sumSWATIperLayer_amt =
        [sum(u_SWATI[soil_discr_to_plot.HorizonNr .== horizon.HorizonNr, :]; dims = [1])'
            for horizon in eachrow(soil_horizons)]
    horizon_labels =
        ["Horizon-$(horizon.HorizonNr): ($(round(horizon.Horizon_Upper_m;digits=2)) to $(round(horizon.Horizon_Lower_m;digits=2)) m)"
        for horizon in eachrow(soil_horizons)]

    # 2) set up the common arguments for all plots below
    # if (compartments == :aboveground || compartments == :belowground)
    #     layout --> RecipesBase.grid(2, 1, heights=[0.2 ,0.8])           #layout --> (2,1)
    #     size --> (1000,700)
    #     idx_d18O_PREC = 1
    #     idx_d18O_SWAT = 2
    #     idx_d2H_PREC  = 1
    #     idx_d2H_SWAT  = 2
    # elseif (compartments == :above_and_belowground)
        # lay = RecipesBase.@layout [RecipesBase.grid(11, 1)]
        lay = RecipesBase.@layout([ °{0.80w} _ ; #1 have no colorbar
                                    °{0.80w} _ ; #2 have no colorbar
                                    °{0.80w} _ ; #3 have no colorbar
                                    °{0.80w} _ ; #4 have no colorbar
                                    °{0.80w} _ ; #5 have no colorbar
                                    ° ; #6
                                    ° ; #7
                                    ° ; #8
                                    ° ;]) #9
        layout --> lay

        size --> (1000,2100)
    #     idx_d18O_PREC = 1
    #     idx_d18O_SWAT = 2
    #     idx_d2H_PREC  = 3
    #     idx_d2H_SWAT  = 4
    # else
    #     # nothing
    # end
    dpi --> 300
    xlim --> xlimits
    leftmargin --> 15mm
    rightmargin --> 15mm

    # # 3) generate plots
    # # NOTE: --> sets attributes only when they don't already exist
    # # NOTE: :=  sets attributes even when they already exist

    # 3a) Precipitation
    @series begin
        # title := "Precipitation"
        seriestype := :bar
        # linecolor := :match
        # line_z := reshape(row_PREC_d18O,1,:)
        # fill_z := reshape(row_PREC_d18O,1,:)
        # clims := clims_d18O
        # colorbar_title := "δ18O [‰]"
        # colorbar := true_to_check_colorbar # TODO: define this once for all plots (except force it for colorbar plot representing common legend)
        yguide := "PREC [mm]"
        legend := false
        subplot := 1
        # and other arguments:
        reshape(x1,1,:), reshape(row_PREC_amt_dense,1,:)
        # reshape(x1_input,1,:), reshape(row_PREC_amt_dense_input,1,:)
    end
    if (compartments == :aboveground || compartments == :above_and_belowground)
        # plot(x, [col_INTS_amt col_INTR_amt col_SNOW_amt col_GWAT_amt col_XYL_amt],
        #      labels = ["INTS" "INTR" "SNOW" "GWAT" "XYL"], ylab = "Amount [mm]")
        # somehow we need a for loop to get the labelling in the correct order
        aboveground_labels = ["INTS" "INTR" "SNOW" "GWAT" "XYL"]
        aboveground_values = [col_INTS_amt col_INTR_amt col_SNOW_amt col_GWAT_amt col_XYL_amt]
        for it in 1:5
            @series begin
                title := "Aboveground"
                ylab := "Amount [mm]"
                subplot := 2
                seriestype := :line
                bg_legend --> colorant"rgba(100%,100%,100%,0.8)"; legend := :topright
                labels := aboveground_labels[it]
                x, aboveground_values[:, it]
            end
        end
    end
    if (compartments == :belowground || compartments == :above_and_belowground)
        # somehow we need a for loop to get the labelling in the correct order
        belowground_labels = ["total SWAT" permutedims(horizon_labels)]
        belowground_values = [col_sumSWATI_amt[:] reduce(hcat, cols_sumSWATIperLayer_amt)]
        for it in 1:length(belowground_labels)
            @series begin
                title := "Belowground"
                ylab := "Amount [mm]"
                subplot := 3
                seriestype := :line
                bg_legend --> colorant"rgba(100%,100%,100%,0.8)"; legend := :topleft
                labels := belowground_labels[it]
                x, belowground_values[:, it]
            end
        end

        # WB error as defined by Ireson 2023
        belowground2_labels = ["Inflow SWATI" "Outflow SWATI" "Storage SWATI" "Inflow Model" "Outflow Model" "Storage Model"]
        belowground2_values = [cum_qIn_SWATI   cum_qOut_SWATI  SWATI_total     cum_qIn_ALL    cum_qOut_ALL    ALL_total]
        belowground2_linestyles = [:solid :solid :solid :dash :dash :dash]
        belowground2_colors = [1 2 3 1 2 3]
        for it in 1:length(belowground2_labels)
            @series begin
                # title := "Belowground"
                # ylab := "Storage [mm] or \n Cumulative In-/Outflow [mm]\n(Soilwater, [+ Groundwater + Temporary Storage])"
                ylab := "Cumulative\nIn-/Outflow or Storage [mm]\nof control volume"
                subplot := 4
                color := belowground2_colors[it]
                linestyle := belowground2_linestyles[it]
                seriestype := :line
                bg_legend --> colorant"rgba(100%,100%,100%,0.8)"; legend := :topleft
                labels := belowground2_labels[it]
                x, belowground2_values[:, it]
            end
        end
        belowground3_labels = ["εB SWATI" "εB Model"]
        belowground3_values = [εB_cumulative_SWATI εB_cumulative_ALL]
        belowground3_linestyles = [:solid :dash]
        for it in 1:length(belowground3_labels)
            @series begin
                # title := "Belowground"
                ylab := "Cumulative\nWater Balance\nError εB [mm]"
                subplot := 5
                color := 3
                linestyle := belowground3_linestyles[it]
                seriestype := :line
                bg_legend --> colorant"rgba(100%,100%,100%,0.8)"; legend := :topleft
                labels := belowground3_labels[it]
                x, belowground3_values[:, it]
            end
        end
        # rows_SWAT_amt0 = u_aux_θ
        rows_SWAT_amt1 = u_SWATI ./solu.prob.p.p_soil.p_THICK  # mm per mm of soil thickness
        rows_SWAT_amt2 = u_SWATI ./ solu.prob.p.p_soil.p_THICK ./ (1 .- solu.prob.p.p_soil.p_STONEF)
        # mm per mm of fine soil thickness (assuming gravel fraction contains no water)
        rows_ψₘpF  = log10.(u_aux_PSIM  .* -10) #(from kPa to log10(hPa))
        rows_ψₜₒₜpF = log10.(u_aux_PSITI .* -10) #(from kPa to log10(hPa))
        # rows_ψₘpF= log10.(u_aux_PSIM  .* -10 + 0.1) #(from kPa to log10(hPa) and small shift to start at pF = -1 instead of pF = -∞)

        cbar_title = "ψₘ [kPa]"
        @series begin
            seriestype := :heatmap
            yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
            colorbar := true_to_check_colorbar; # clims := clims_d2H
            yguide := "Depth [mm]"; colorbar_title := cbar_title
            c := cgrad_two_separate_gradients_at_mid(u_aux_PSIM, -0.00001) #shows 0 in red
            # c := cgrad_two_separate_gradients_at_mid(u_aux_PSIM, -0.0)     # shows 0 in blue
            subplot := 6
            # and other arguments:
            x, y_center, u_aux_PSIM;
        end
        # cbar_title = "ψₜ [kPa]"
        # @series begin
        #     seriestype := :heatmap
        #     yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
        #     colorbar := true_to_check_colorbar; # clims := clims_d2H
        #     yguide := "Depth [mm]"; colorbar_title := cbar_title
        #     c := cgrad(color_scheme,  rev = true)
        #     subplot := 6
        #     # and other arguments:
        #     x, y_center, u_aux_PSITI;
        # end
        # cbar_title = log₁₀(-ψₘ hPa)"
        # @series begin
        #     seriestype := :heatmap
        #     yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
        #     colorbar := true_to_check_colorbar; # clims := clims_d2H
        #     yguide := "Depth [mm]"; colorbar_title := cbar_title
        #     # c := cgrad(color_scheme,  rev = true)
        #     c := cgrad_two_separate_gradients_at_mid(rows_ψₘpF,)
        #     subplot := 6
        #     # and other arguments:
        #     x, y_center, rows_ψₘpF
        # end
        if (RWUcentroid == :showRWUcentroid)
            @series begin
                #plot!(x, row_RWU_centroid_mm', yflip=true, color=:white, label = "")
                color := :white
                label := RWUcentroidLabel
                #bg_legend --> colorant"rgba(100%,100%,100%,0.8)";legend := :bottomright
                bg_legend --> colorant"rgba(100%,100%,100%,0.0)"; legend := :bottomright; fg_legend --> :transparent; legendfontcolor := :white
                yflip := true; yticks := y_soil_ticks
                yguide := "Depth [mm]"; colorbar_title := cbar_title
                subplot := 6
                x, row_RWU_centroid_mm'
            end
        end
        # @series begin
        #     seriestype := :heatmap
        #     yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
        #     colorbar := true_to_check_colorbar; # clims := clims_d2H
        #     yguide := "Depth [mm]"; colorbar_title := "pF = log₁₀(-ψₜₒₜ hPa)" #colorbar_title := "pF = \nlog₁₀(-ψ hPa)"
        #     c := cgrad(color_scheme,  rev = true)
        #     subplot := 7
        #     # and other arguments:
        #     x, y_center, rows_ψₜₒₜpF
        # end
        # @series begin
        #     seriestype := :heatmap
        #     yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
        #     colorbar := true_to_check_colorbar; # clims := clims_d2H
        #     yguide := "Depth [mm]"; colorbar_title := "SWATI [mm]"
        #     c := cgrad(color_scheme,  rev = false)
        #     subplot := 7
        #     # and other arguments:
        #     x, y_center, u_SWATI;                 # deactivated u_SWATI as it is resolution dependent!
        # end
        @series begin
            seriestype := :heatmap
            yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
            colorbar := true_to_check_colorbar; # clims := clims_d2H
            yguide := "Depth [mm]"; colorbar_title := "Wetness [-]"
            # c := cgrad(color_scheme,  rev = false)
            # c := cgrad_two_separate_gradients_at_mid(u_aux_WETNES, 1.0; down = 0.0, up = 1.0); clims := (0,1)
            # c := cgrad_two_separate_gradients_at_mid(u_aux_WETNES, 1.0; down = 0.0, up = 1.0)
            c := cgrad_two_separate_gradients_at_mid(u_aux_WETNES, 1.0)
            subplot := 7
            # and other arguments:
            x, y_center, u_aux_WETNES;
        end
        @series begin
            seriestype := :heatmap
            yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
            colorbar := true_to_check_colorbar; # clims := clims_d2H
            yguide := "Depth [mm]"; colorbar_title := "θ [m3/m3]\n(of fine soil volume)" #colorbar_title := "θ [m3/m3]\n(fine soil)" # "θ [m3/m3]"
            c := cgrad(color_scheme,  rev = false)
            subplot := 8
            # and other arguments:
            x, y_center, u_aux_θ;
        end
        # rows_RWU = rows_RWU_mmDay ./ solu.prob.p.p_soil.p_THICK
        # @series begin
        #     seriestype := :heatmap
        #     yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
        #     colorbar := true_to_check_colorbar; # clims := clims_d2H
        #     yguide := "Depth [mm]"; colorbar_title := "RWU [mm water/day per mm soil depth]"
        #     c := :diverging_bwr_20_95_c54_n256; clim := maximum(abs.(rows_RWU)) .* (-1, 1)
        #     subplot := 8
        #     # and other arguments:
        #     x, y_center, rows_RWU;
        # end
        @series begin
            seriestype := :heatmap
            yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
            colorbar := true_to_check_colorbar; # clims := clims_d2H
            yguide := "Depth [mm]"; colorbar_title := "θ [m3/m3]\n(of total volume, incl stonef)" #colorbar_title := "θ [-]\n(total, incl stonef)"
            # c := cgrad(color_scheme,  rev = false)
            # c := cgrad_two_separate_gradients_at_mid(rows_SWAT_amt1, 1.0; down = 0.0, up = 1.0); clims := (0,1)
            # c := cgrad_two_separate_gradients_at_mid(rows_SWAT_amt1, 1.0; down = 0.0, up = 1.0)
            c := cgrad_two_separate_gradients_at_mid(rows_SWAT_amt1, 1.0)
            subplot := 9
            # and other arguments:
            x, y_center, rows_SWAT_amt1;
        end
        # @series begin
        #     seriestype := :heatmap
        #     yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
        #     colorbar := true_to_check_colorbar
        #     yguide := "Depth [mm]"; colorbar_title := "θ [-] (fine soil 2)"#colorbar_title := "θ [-]\n(fine soil 2)"
        #     c := cgrad(color_scheme,  rev = false)
        #     subplot := 10
        #     # and other arguments:
        #     x, y_center, rows_SWAT_amt2           # deactivated: as it was same as `x, y_center, u_aux_θ;`
        # end
        # @series begin
        #     seriestype := :heatmap
        #     yflip := true; yticks := y_soil_ticks #(y_ticks, y_labels)
        #     colorbar := true_to_check_colorbar; # clims := clims_d2H
        #     yguide := "Depth [mm]"; colorbar_title := "K [mm/day]"
        #     c := cgrad(color_scheme,  rev = false)
        #     subplot := 10
        #     # and other arguments:
        #     x, y_center, p_fu_KK;
        #     # x, y_center, log10.(p_fu_KK);
        # end

        # TODO: edges of cells in heatmap are not entirely correct. Find a way to override heatmap()
        #       where we provide cell edges (n+1) instead of cell centers (n)
        #       e.g. plots_heatmap_edges: @recipe function f(::Type{Val{:plots_heatmap_edges}}, xe, ye, z)
        #       e.g. plots_heatmap_edges:     m, n = size(z.surf)
        #       e.g. plots_heatmap_edges:     x_pts, y_pts = fill(NaN, 6 * m * n), fill(NaN, 6 * m * n)
        #       e.g. plots_heatmap_edges:     fz = zeros(m * n)
        #       e.g. plots_heatmap_edges:     for i in 1:m # y
        #       e.g. plots_heatmap_edges:         for j in 1:n # x
        #       e.g. plots_heatmap_edges:             k = (j - 1) * m + i
        #       e.g. plots_heatmap_edges:             inds = (6 * (k - 1) + 1):(6 * k - 1)
        #       e.g. plots_heatmap_edges:             x_pts[inds] .= [xe[j], xe[j + 1], xe[j + 1], xe[j], xe[j]]
        #       e.g. plots_heatmap_edges:             y_pts[inds] .= [ye[i], ye[i], ye[i + 1], ye[i + 1], ye[i]]
        #       e.g. plots_heatmap_edges:             fz[k] = z.surf[i, j]
        #       e.g. plots_heatmap_edges:         end
        #       e.g. plots_heatmap_edges:     end
        #       e.g. plots_heatmap_edges:     ensure_gradient!(plotattributes, :fillcolor, :fillalpha)
        #       e.g. plots_heatmap_edges:     fill_z := fz
        #       e.g. plots_heatmap_edges:     line_z := fz
        #       e.g. plots_heatmap_edges:     x := x_pts
        #       e.g. plots_heatmap_edges:     y := y_pts
        #       e.g. plots_heatmap_edges:     z := nothing
        #       e.g. plots_heatmap_edges:     seriestype := :shape
        #       e.g. plots_heatmap_edges:     label := ""
        #       e.g. plots_heatmap_edges:     widen --> false
        #       e.g. plots_heatmap_edges:     ()
        #       e.g. plots_heatmap_edges: end
        #       e.g. plots_heatmap_edges: @deps plots_heatmap_edges shape
        #       e.g. plots_heatmap_edges: @shorthands plots_heatmap_edges
        #       e.g. plots_heatmap_edges:
        #       e.g. plots_heatmap_edges: Plots.heatmap(x[1:100], y_center, z[:,1:100])
        #       e.g. plots_heatmap_edges: Plots.heatmap(x[1:100], y_center, z[:,1:100])
        #       e.g. plots_heatmap_edges: plot(t = :heatmap, x[1:50], y_center, z[:,1:50]) # works
        #       e.g. plots_heatmap_edges: plot(t = :plots_heatmap, x[1:50], y_center, z[:,1:50]) # doesn't work
        #       e.g. plots_heatmap_edges: plot(t = :plots_heatmap_edges, x[1:50], y_center, z[:,1:50]) # doesn't work either
    end
end

"""
    plotisotopes(simulation::DiscretizedSPAC)
    plotisotopes(simulation::DiscretizedSPAC, isotope::Symbol)
    plotisotopes(simulation::DiscretizedSPAC, isotope::Symbol, (d18O = (-6, -16), d2H = (-125, -40)))
    plotisotopes(simulation::DiscretizedSPAC, isotope::Symbol, (d18O = (-6, -16), d2H = (-125, -40)), RWUcentroid::Symbol))

Plots the isotope results of a SPAC Simulation. By default both δ18O and δ2H.
The user can override this with the second argument isotope as one of `:d18O`, `:d2H`, or `:d18O_and_d2H`.
RWUcentroid can have values of either `:dontShowRWUcentroid` or `:showRWUcentroid`.
"""
@userplot PlotIsotopes
@recipe function f(pliso::PlotIsotopes)
    # 0) parse input arguments
    if length(pliso.args) == 4
        simulation = pliso.args[1]
        isotope    = pliso.args[2]
        clims      = pliso.args[3]
        RWUcentroid= pliso.args[4]
    elseif length(pliso.args) == 3
        simulation = pliso.args[1]
        isotope    = pliso.args[2]
        clims      = pliso.args[3]
        RWUcentroid= :dontShowRWUcentroid
    elseif length(pliso.args) == 2
        simulation = pliso.args[1]
        isotope    = pliso.args[2]
        clims      = (d18O = :auto, d2H = :auto)
        RWUcentroid= :dontShowRWUcentroid
    elseif length(pliso.args) == 1
        simulation = pliso.args[1]
        isotope    = :d18O_and_d2H
        clims      = (d18O = :auto, d2H = :auto)
        RWUcentroid= :dontShowRWUcentroid
    else
        error("plotisotopes requires an unnamed first argument of type DiscretizedSPAC, and optional unnamed second/third arguments (:d18O_and_d2H, :d18O, or :d2H) and (:dontShowRWUcentroid, :showRWUcentroid). Other arguments to plot() should be separated by `;`.")
    end
    if !(RWUcentroid == :dontShowRWUcentroid || RWUcentroid == :showRWUcentroid)
        error("Third unnamed argument to plotisotopes should be one of (:dontShowRWUcentroid, :showRWUcentroid). Got: $(RWUcentroid)")
    end
    if !(isotope == :d18O || isotope == :d2H || isotope == :d18O_and_d2H)
        error("Second unnamed argument to plotisotopes should be one of (:d18O_and_d2H, :d18O, or :d2H). Got: $(isotope)")
    end
    if !(simulation isa DiscretizedSPAC)
        error("First unnamed argument to plotisotopes should be of type DiscretizedSPAC. Got: $(typeof(simulation))")
    end
    if isnothing(simulation.ODESolution)
        error("plotisotopes requires a solved system. Please `simulate!()` the DiscretizedSPAC first.")
    end

    # 1) prepare data to plot
    solu = simulation.ODESolution

    # 1a) extract data from solution object `solu`
    df, RWUcentroidLabel = intern___get_data_for_isotopePlot(simulation)

    days_to_read_out_d = df.days_to_read_out_d
    t_ref = solu.prob.p.REFERENCE_DATE
    x = RelativeDaysFloat2DateTime.(days_to_read_out_d, t_ref);
    y_center = cumsum(solu.prob.p.p_soil.p_THICK) - solu.prob.p.p_soil.p_THICK/2

    simulate_isotopes = solu.prob.p.simulate_isotopes
    @assert simulate_isotopes "Provided DiscretizedSPAC() did not simulate isotopes"

    simulate_irrigation = solu.prob.p.simulate_irrigation
    
    # Some hardcoded options:
    xlimits = RelativeDaysFloat2DateTime.(solu.prob.tspan, t_ref)
    tick_function = (x1, x2) -> PlotUtils.optimize_ticks(x1, x2; k_min = 4)
    # color_scheme = :default # https://docs.juliaplots.org/latest/generated/colorschemes
    # color_scheme = :blues
    color_scheme = :heat

    # Results
    maxdepth = maximum(cumsum(solu.prob.p.p_soil.p_THICK))
    t_ref = solu.prob.p.REFERENCE_DATE
    t1 = range(extrema(solu.t)..., step = 1) # Plot forcing as daily, even if solution output (ODESolution.t) is not dense
    x1 = RelativeDaysFloat2DateTime.(t1, t_ref);

    # Make rows for plotting as heatmap
    row_PREC_amt_dense  = reshape(df.col_PREC_amt_dense, 1, :)
    row_PREC_d18O_dense = reshape(df.col_PREC_d18O_dense, 1, :)
    row_PREC_d2H_dense  = reshape(df.col_PREC_d2H_dense, 1, :)

    rows_SWAT_d18O = Matrix(permutedims(df[:, r"d18O_[0-9]"] ))
    rows_SWAT_d2H  = Matrix(permutedims(df[:, r"d2H_[0-9]"] ))
    row_NaN             = fill(NaN, 1,size(df,1))
    row_PREC_d18O = reshape(df.col_PREC_d18O, 1, :)
    row_PREC_d2H  = reshape(df.col_PREC_d2H, 1, :)
    row_INTS_d18O = reduce(hcat, df.col_INTS_d18O)
    row_INTR_d18O = reduce(hcat, df.col_INTR_d18O)
    row_SNOW_d18O = reduce(hcat, df.col_SNOW_d18O)
    row_IRRIG_d18O = reduce(hcat, df.col_IRRIG_d18O)
    row_GWAT_d18O = reduce(hcat, df.col_GWAT_d18O)
    row_RWU_d18O = reduce(hcat, df.col_RWU_d18O)
    row_XYL_d18O = reduce(hcat, df.col_XYL_d18O)
    row_INTS_d2H = reduce(hcat, df.col_INTS_d2H)
    row_INTR_d2H = reduce(hcat, df.col_INTR_d2H)
    row_SNOW_d2H = reduce(hcat, df.col_SNOW_d2H)
    row_IRRIG_d2H = reduce(hcat, df.col_IRRIG_d2H)
    row_GWAT_d2H = reduce(hcat, df.col_GWAT_d2H)
    row_RWU_d2H = reduce(hcat, df.col_RWU_d2H)
    row_XYL_d2H = reduce(hcat, df.col_XYL_d2H)

    row_RWU_centroid_mm = permutedims(df.col_RWU_centroid_mm)

    # 1b) define some plot arguments based on the extracted data
    # color scheme:
    true_to_check_colorbar = false; # set this flag to false for final plot, true for debugging.

    clims_d18O = ifelse(clims.d18O == :auto, extrema(filter(!isnan, row_PREC_d18O)), clims.d18O)
    clims_d2H  = ifelse(clims.d2H  == :auto, extrema(filter(!isnan, row_PREC_d2H)),  clims.d2H)
    @show clims_d18O
    @show clims_d2H
    # 2) set up the common arguments for all plots below
    if (isotope == :d18O || isotope == :d2H)
        lay = RecipesBase.@layout([ °{0.91w} _ ;   #1 have no colorbar
                                    °          ;]) #2
        layout --> lay
        size --> (1000,700)
        idx_d18O_PREC = 1
        idx_d18O_SWAT = 2
        idx_d2H_PREC  = 1
        idx_d2H_SWAT  = 2
    elseif (isotope == :d18O_and_d2H)
        lay = RecipesBase.@layout([ °{0.86w} _ ;   #1 have no colorbar
                                    °          ;   #2
                                    °{0.86w} _ ;   #3 have no colorbar
                                    °          ;]) #4
        layout --> lay
        size --> (1000,1400)
        idx_d18O_PREC = 1
        idx_d18O_SWAT = 2
        idx_d2H_PREC  = 3
        idx_d2H_SWAT  = 4
    else
        # nothing
    end
    dpi --> 300
    xlim --> xlimits
    leftmargin --> 15mm

    # # 3) generate plots
    # # NOTE: --> sets attributes only when they don't already exist
    # # NOTE: :=  sets attributes even when they already exist

    # 3a) Precipitation
    # We reproduce the following plot:
    # Workaround for barplot # https://github.com/JuliaPlots/Plots.jl/issues/3880
    # ts_PREC_δ18O = plot(reshape(x,1,:), reshape(row_PREC_amt_dense,1,:), t = :bar, line_z = reshape(row_PREC_d18O,1,:), fill_z = reshape(row_PREC_d18O,1,:),
    #                 clims=clims_d18O, colorbar = true_to_check_colorbar, legend = false,
    #                 ylabel = "PREC [mm]"); # TODO(bernhard): make this work with a barplot
    # ts_PREC_δ2H = ...
    if (isotope == :d18O || isotope == :d18O_and_d2H)
        @series begin # ts_PREC_δ18O
            title --> "δ18O"
            seriestype := :bar
            # linecolor := :match
            line_z := reshape(row_PREC_d18O_dense,1,:)
            fill_z := reshape(row_PREC_d18O_dense,1,:)
            clims := clims_d18O
            colorbar_title := "δ18O [‰]"
            yguide := "PREC [mm]"
            colorbar := true_to_check_colorbar #https://stackoverflow.com/a/59257011
            legend := false;
            subplot := idx_d18O_PREC

            # and other arguments:
            reshape(x1,1,:), reshape(row_PREC_amt_dense,1,:)
        end
    end
    if (isotope == :d2H || isotope == :d18O_and_d2H)
        @series begin # ts_PREC_δ2H
            title --> "δ2H"
            seriestype := :bar
            # linecolor := :match
            line_z := reshape(row_PREC_d2H_dense,1,:)
            fill_z := reshape(row_PREC_d2H_dense,1,:)
            clims := clims_d2H
            colorbar_title := "δ2H [‰]"
            yguide := "PREC [mm]"
            colorbar := true_to_check_colorbar #https://stackoverflow.com/a/59257011
            legend := false;
            subplot := idx_d2H_PREC

            # and other arguments:
            reshape(x1,1,:), reshape(row_PREC_amt_dense,1,:)
        end
    end

    # 3b) Heatmap (containing SWATI and other compartments)
    # y_labels   = ["INTS"; ""; "INTR"; ""; "SNOW"; ""; round.(y_center); "";             "GWAT"]
    # y_soil_ticks = optimize_ticks(extrema(y_center)...; k_min = 4)[1]
    # y_soil_ticks = optimize_ticks(0., round(maxdepth))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
    y_extended = [-500; -350; -300; -250; -200; -150; -100; -50;         y_center;             (maxdepth .+ [50; 100; 150; 250; 300;400])]
    y_soil_ticks = tick_function(0., round(maxdepth))[1] # TODO(bernhard): how to do without loading Plots.optimize_ticks()
    y_ticks    = [-500;       -300;       -200;       -100;          y_soil_ticks;             (maxdepth .+ [    100;      250;     400])]
    
    if simulate_irrigation # replace snow isotopic signature with irrigation input
        y_labels   = ["PREC";   "INTS";     "INTR";     "IRRIG";     round.(y_soil_ticks; digits=0);                                                "GWAT";    "RWU";     "XYLEM"]
        z2_extended = [row_PREC_d18O; row_NaN; row_INTS_d18O; row_NaN; row_INTR_d18O; row_NaN; row_IRRIG_d18O; row_NaN; rows_SWAT_d18O; row_NaN; row_GWAT_d18O; row_NaN; row_RWU_d18O; row_NaN; row_XYL_d18O]
        z3_extended = [row_PREC_d2H;  row_NaN; row_INTS_d2H;  row_NaN; row_INTR_d2H;  row_NaN; row_IRRIG_d2H;  row_NaN; rows_SWAT_d2H;  row_NaN; row_GWAT_d2H;  row_NaN; row_RWU_d2H;  row_NaN; row_XYL_d2H ]
    else
        y_labels   = ["PREC";   "INTS";     "INTR";     "SNOW";     round.(y_soil_ticks; digits=0);                                                "GWAT";    "RWU";     "XYLEM"]
        z2_extended = [row_PREC_d18O; row_NaN; row_INTS_d18O; row_NaN; row_INTR_d18O; row_NaN; row_SNOW_d18O; row_NaN; rows_SWAT_d18O; row_NaN; row_GWAT_d18O; row_NaN; row_RWU_d18O; row_NaN; row_XYL_d18O]
        z3_extended = [row_PREC_d2H;  row_NaN; row_INTS_d2H;  row_NaN; row_INTR_d2H;  row_NaN; row_SNOW_d2H;  row_NaN; rows_SWAT_d2H;  row_NaN; row_GWAT_d2H;  row_NaN; row_RWU_d2H;  row_NaN; row_XYL_d2H ]
    end

    if (isotope == :d18O || isotope == :d18O_and_d2H)
        @series begin # pl_δ18O
            seriestype := :heatmap
            yflip := true
            yticks := (y_ticks, y_labels)
            yguide := "Depth [mm]"
            colorbar_title := "δ18O [‰]"
            clims := clims_d18O
            colorbar := true # colorbar := true_to_check_colorbar
            subplot := idx_d18O_SWAT

            # and other arguments:
            x, y_extended, z2_extended;
        end
        if (RWUcentroid == :showRWUcentroid)
            @series begin
                color := :white
                label := "mean RWU depth"
                bg_legend --> colorant"rgba(100%,100%,100%,0.0)"; legend := :right;       fg_legend --> :transparent; legendfontcolor := :white
                # bg_legend --> colorant"rgba(100%,100%,100%,0.0)"; legend := :bottomright; fg_legend --> :transparent; legendfontcolor := :white
                yflip := true; yticks := (y_ticks, y_labels)
                yguide := "Depth [mm]"; colorbar_title := "δ18O [‰]"
                subplot := idx_d18O_SWAT
                x, row_RWU_centroid_mm'
            end
        end
    end

    if (isotope == :d2H || isotope == :d18O_and_d2H)
        @series begin # pl_δ2H
            seriestype := :heatmap
            yflip := true
            yticks := (y_ticks, y_labels)
            yguide := "Depth [mm]"
            colorbar_title := "δ2H [‰]"
            clims := clims_d2H
            colorbar := true # colorbar := true_to_check_colorbar
            subplot := idx_d2H_SWAT

            # and other arguments:
            x, y_extended, z3_extended;
        end
        if (RWUcentroid == :showRWUcentroid)
            @series begin
                color := :white
                label := RWUcentroidLabel
                bg_legend --> colorant"rgba(100%,100%,100%,0.0)";legend := :right; fg_legend --> :transparent; legendfontcolor := :white
                #bg_legend --> colorant"rgba(100%,100%,100%,0.0)"; legend := :right; fg_legend --> :transparent; legendfontcolor := :white
                yflip := true; yticks := (y_ticks, y_labels)
                yguide := "Depth [mm]"; colorbar_title := "δ2H [‰]"
                subplot := idx_d2H_SWAT
                x, row_RWU_centroid_mm'
            end
        end
    end

    # 3c) Colorbar
    # ....
    # TODO: this code needs to reproduce below steps:
    # pl_colorbar_δ18O = plot([0,0], [0,1], zcolor=[0,1], t=:scatter, xlims=(1,1.1), # set xlims that no marker is shown
    #                         clims=clims_d18O, colorbar_title="δ180 [‰]",
    #                         grid=false, showaxis=false, ticks=false, label=false);
    # pl_colorbar_δ2H  = plot([0,0], [0,1], zcolor=[0,1], t=:scatter, xlims=(1,1.1), # set xlims that no marker is shown
    #                         clims=clims_d18O, colorbar_title="δ2H [‰]",
    #                         grid=false, showaxis=false, ticks=false, label=false);

    # l = @layout [grid(2, 1, heights=[0.2, 0.8]) a{0.055w}]
    # pl_final_δ18O = plot(ts_PREC_δ18O, pl_δ18O, pl_colorbar_δ18O,
    #     clims = clims_d18O,
    #     layout = l, link = :x);
    # l = @layout [grid(2, 1, heights=[0.2, 0.8]) a{0.055w}]
    # pl_final_δ2H = plot(ts_PREC_δ2H,  pl_δ2H, pl_colorbar_δ2H,
    #     clims = clims_d2H,
    #     layout = l, link = :x);
end



"""
    plotforcingandstates(simulation::DiscretizedSPAC)

Plots the forcing, states and major fluxes as results of a SPAC Simulation.
"""
@userplot PlotForcingAndStates
@recipe function f(plfor::PlotForcingAndStates)
    # 0) parse input arguments
    if length(plfor.args) == 1
        simulation = plfor.args[1]
    else
        error("plotforcingandstates requires an unnamed first argument of type DiscretizedSPAC. Other arguments to plot() should be separated by `;`.")
    end
    if !(simulation isa DiscretizedSPAC)
        error("First unnamed argument to plotamounts should be of type DiscretizedSPAC. Got: $(typeof(simulation))")
    end
    if isnothing(simulation.ODESolution)
        error("plotforcingandstates requires a solved system. Please `simulate!()` the DiscretizedSPAC first.")
    end

    # 1) prepare data to plot
    solu = simulation.ODESolution;

    # 1a) extract data from solution object `solu`
    # 1a) i) forcing (wind, vappres, globrad, tmax, tmin, prec, lai, )
    t_ref = solu.prob.p.REFERENCE_DATE
    # t1 = range(simulation.ODEProblem.tspan..., step = 1)           # Plot forcing as daily, even if solution output (ODESolution.t) is not dense
    # t1 = range(extrema(simulation.ODESolution.t)..., step = 1)     # Plot forcing as daily, even if solution output (ODESolution.t) is not dense
    t1 = simulation.parametrizedSPAC.forcing.meteo["p_days"]         # forcing period can be longer than simulation output (e.g. when using a spinup period)
    x1 = RelativeDaysFloat2DateTime.(t1, t_ref);
    y11 = simulation.ODESolution.prob.p.p_WIND.(t1);      lbl11 = "p_WIND [m/s]";
    y12 = simulation.ODESolution.prob.p.p_VAPPRES.(t1);   lbl12 = "p_VAPPRES [kPa]";
    y13 = simulation.ODESolution.prob.p.p_GLOBRAD.(t1);   lbl13 = "p_GLOBRAD [MJ/day/m2]"
    y14 = hcat(simulation.ODESolution.prob.p.p_TMIN.(t1),
               simulation.ODESolution.prob.p.p_TMAX.(t1));lbl14 = ["p_TMIN [°C]" "p_TMAX [°C]"]
    y15 = simulation.ODESolution.prob.p.p_PREC.(t1);      lbl15 = "p_PREC [mm]"
    # y15 = simulation.ODESolution.prob.p.p_IRRIG.(t1);      lbl15 = "p_IRRIG [mm]"
    # plot_forcing = plot(layout = (:,1),
    #     plot(x1, y11; labels = lbl11),
    #     plot(x1, y12; labels = lbl12),
    #     plot(x1, y13; labels = lbl13),
    #     plot(x1, y14; labels = lbl14),
    #     plot(x1, y15; labels = lbl15))

    t2 = t1; x2 = x1;
    y21 = hcat(simulation.ODESolution.prob.p.p_DENSEF.(t2),
                        simulation.ODESolution.prob.p.p_SAI.(t2),
                        simulation.ODESolution.prob.p.p_LAI.(t2)); lbl21 = ["p_DENSEF [-]" "p_SAI [-]" "p_LAI [-]"];
    y22 = simulation.ODESolution.prob.p.p_HEIGHT.(t2);             lbl22 = "p_HEIGHT [-]";
    # plot_vegetation = plot(layout = (2,1), title = "Vegetation",
    #     plot(x2, y21; labels = lbl21),
    #     plot(x2, y22;labels = lbl22))
    #     # plot!(twinx(),100*rand(10))
    #     # plot!(twinx(), simulation.ODESolution_datetime, simulation.ODESolution.prob.p.p_HEIGHT.(t2))               #;labels = "p_HEIGHT [-]")

    # 1a) ii) states (scalar and vector)
    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) = LWFBrook90.intern___get_SWATI_derivatives(simulation)

    x3 = simulation.ODESolution_datetime;
    y31 = hcat(reduce(hcat, [simulation.ODESolution(t).INTR.mm   for t in simulation.ODESolution.t])',
                reduce(hcat, [simulation.ODESolution(t).INTS.mm   for t in simulation.ODESolution.t])',
                reduce(hcat, [simulation.ODESolution(t).SNOW.mm   for t in simulation.ODESolution.t])',
                reduce(hcat, [simulation.ODESolution(t).SNOWLQ.mm   for t in simulation.ODESolution.t])',
                reduce(hcat, [simulation.ODESolution(t).XYLEM.mm   for t in simulation.ODESolution.t])');
    lbl31 = ["INTR [mm]" "INTS [mm]" "SNOW [mm]" "SNOWLQ [mm]" "XYLEM [mm]"];
    y32 = hcat(sum(u_SWATI; dims=1)',
                reduce(hcat, [simulation.ODESolution(t).GWAT.mm   for t in simulation.ODESolution.t])');
    lbl32 = ["Total Soil Water [mm]" "GWAT [mm]"];
    # pl_states1 = plot(x3,
    #         # yaxis=:log,  1 .+
    #         y31; title = "Scalar states", labels = lbl31)
    # pl_states2 = plot(x3, y32 = ; title = "Soil/ground water", labels = lbl32)
    # plot_states = plot(pl_states1, pl_states2, layout = (2,1))

    # 1a) iii) fluxes
    x4 = simulation.ODESolution_datetime;
    y41 = hcat(  reduce(hcat, [simulation.ODESolution(t).RWU.mmday   for t in simulation.ODESolution.t])',
                                        sum(reduce(hcat, [simulation.ODESolution(t).TRANI.mmday   for t in simulation.ODESolution.t]); dims=1)')
    lbl41 = ["RWU (mm/day)" "sum(TRANI)"]
    # plot_fluxes = plot(x4, y41; labels = lbl41)

    # plot(plot_forcing, plot_vegetation, plot_states, plot_fluxes,
    #     layout = grid(4, 1, heights=[0.5 ,0.2, 0.2, 0.1]),
    #     size=(800,2000), dpi = 300, leftmargin = 15mm)

    # 2) set up the common arguments for all plots below
    # lay = RecipesBase.grid(4, 1, heights=[0.5 ,0.2, 0.2, 0.1])
    lay = RecipesBase.grid(10, 1)
    layout --> lay
    size --> (800,2000)
    dpi --> 300
    # xlim --> xlimits
    leftmargin --> 15mm
    legend --> :topleft

    # # 3) generate plots
    # # NOTE: --> sets attributes only when they don't already exist
    # # NOTE: :=  sets attributes even when they already exist

    # WORKAROUND:
    # somehow we need a for loop to get the labelling in the correct order
    function own_plotseries(yguide, subplot, labels, x, y)
        for it in 1:size(y,2)
            @series begin
                yguide := yguide
                subplot := subplot
                labels := labels[it]
                x, y[:,it]
            end
        end
    end
    own_plotseries("", 1, [lbl11], x1, y11)
    own_plotseries("", 2, [lbl12], x1, y12)
    own_plotseries("", 3, [lbl13], x1, y13)
    own_plotseries("", 4, lbl14, x1, y14)
    own_plotseries("", 5, [lbl15], x1, y15)
    own_plotseries("", 6, lbl21, x2, y21)
    own_plotseries("", 7, [lbl22], x2, y22)
    own_plotseries("", 8, lbl31, x3, y31)
    own_plotseries("", 9, lbl32, x3, y32)
    own_plotseries("", 10, lbl41, x4, y41)
end
############################################################################################
############################################################################################
############################################################################################