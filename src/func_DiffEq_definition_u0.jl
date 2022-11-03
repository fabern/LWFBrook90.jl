"""
    define_LWFB90_u0()

Generate vector u0 needed for ODE() problem in DiffEq.jl package.
"""
function define_LWFB90_u0(;simulate_isotopes, compute_intermediate_quantities, NLAYER)
    # u0_treespecies_names = (default = 1) # e.g. (beech = 1, oak = 2, spruce = 3)
    names_accum = compute_intermediate_quantities ? reshape([
                "cum_d_prec";"cum_d_rfal";"cum_d_sfal";"cum_d_rint"; "cum_d_sint";"cum_d_rsno";
                "cum_d_rnet";"cum_d_smlt";"cum_d_evap";"cum_d_tran";"cum_d_irvp";"cum_d_isvp";
                "cum_d_slvp";"cum_d_snvp";"cum_d_pint";"cum_d_ptran";"cum_d_pslvp";
                "flow";"seep";"srfl";"slfl";"byfl";"dsfl";"gwfl";"vrfln";
                "cum_d_rthr";"cum_d_sthr";
                "totalSWAT"; "new_totalWATER"; "BALERD_SWAT"; "BALERD_total";
            ],1,:) : []

    variable_names = simulate_isotopes ? (d18O = 2, d2H = 3) : ()
    N_isotopes             = length(variable_names)
    N_separate_treespecies = 1 # TODO(bernhard) currently only one species is implemented
    N_accum_var            = length(names_accum)

    u_totalRWUinit_mmday = zeros(1,  1+N_isotopes, N_separate_treespecies)
    u_Xyleminit_mm       = zeros(1,  1+N_isotopes, N_separate_treespecies)  #[5.0 -12 -95] .* ones(1, 1+N_isotopes, N_separate_treespecies) # # start out with same concentration as in first soil layer
    u_TRANIinit_mmday    = zeros(NLAYER, 1+N_isotopes, N_separate_treespecies)

    u0_NamedTuple = (GWAT  = zeros(1, 1+N_isotopes, 1),
            INTS   = zeros(1, 1+N_isotopes, 1),
            INTR   = zeros(1, 1+N_isotopes, 1),
            SNOW   = zeros(1, 1+N_isotopes, 1),
            CC     = zeros(1, 1, 1),
            SNOWLQ = zeros(1, 1, 1),
            SWATI  = zeros(NLAYER, 1+N_isotopes, 1), #SWATI  = zeros(p[1][2][1], 1+N_isotopes)) # p[1][2][1] = NLAYER
            RWU    = u_totalRWUinit_mmday,
            XYLEM  = u_Xyleminit_mm,
            TRANI  = u_TRANIinit_mmday,
            # Further structures for auxiliary soil variables (θ,ψ,K) and accumulation variables
            aux    = zeros(NLAYER, 3), # TODO: where to store θ, ψ and K(θ) ?
            accum  = zeros(N_accum_var,1))
    # Give ComponentArray as u0 to DiffEq.jl
    name_variables = ifelse(simulate_isotopes, (:mm, :d18O, :d2H), (:mm,))
    name_aux       = (:θ,:ψ,:K)
    name_accum     = (:cum_d_prec, :cum_d_rfal, :cum_d_sfal, :cum_d_rint,  :cum_d_sint, :cum_d_rsno,
                    :cum_d_rnet, :cum_d_smlt, :cum_d_evap, :cum_d_tran, :cum_d_irvp, :cum_d_isvp,
                    :cum_d_slvp, :cum_d_snvp, :cum_d_pint, :cum_d_ptran, :cum_d_pslvp,
                    :flow, :seep, :srfl, :slfl, :byfl, :dsfl, :gwfl, :vrfln,
                    :cum_d_rthr, :cum_d_sthr,
                    :totalSWAT,  :new_totalWATER,  :BALERD_SWAT,  :BALERD_total)
    if simulate_isotopes
        u0 = ComponentArray(
            SWATI  = NamedTuple{name_variables, NTuple{3, Vector{Float64}}}(tuple(eachcol(u0_NamedTuple[:SWATI][:,:,1])...)),
            GWAT   = NamedTuple{name_variables,      NTuple{3, Float64}}(u0_NamedTuple[:GWAT][:,:,1]),
            INTS   = NamedTuple{name_variables,      NTuple{3, Float64}}(u0_NamedTuple[:INTS][:,:,1]),
            INTR   = NamedTuple{name_variables,      NTuple{3, Float64}}(u0_NamedTuple[:INTR][:,:,1]),
            SNOW   = NamedTuple{name_variables,      NTuple{3, Float64}}(u0_NamedTuple[:SNOW][:,:,1]),
            RWU    = NamedTuple{name_variables,      NTuple{3, Float64}}(u0_NamedTuple[:RWU]        ),
            XYLEM  = NamedTuple{name_variables,      NTuple{3, Float64}}(u0_NamedTuple[:XYLEM]      ),
            CC     = NamedTuple{name_variables[[1]], NTuple{1, Float64}}(u0_NamedTuple[:CC][:,:,1]),
            SNOWLQ = NamedTuple{name_variables[[1]], NTuple{1, Float64}}(u0_NamedTuple[:SNOWLQ][:,:,1]),

            TRANI = NamedTuple{name_variables, NTuple{3, Vector{Float64}}}(tuple(eachcol(u0_NamedTuple[:TRANI][:,:,1])...)),

            aux   = NamedTuple{name_aux, NTuple{3, Vector{Float64}}}(tuple(eachcol(u0_NamedTuple[:aux])...)),
            accum = NamedTuple{name_accum, NTuple{31, Float64}}((0. for i in eachindex(name_accum))))
    else
        # TODO(bernhard): check if this is bad programming if NTuple{1, ...} depends on runtime variable simulate_isotopes...
        u0 = ComponentArray(
            SWATI  = NamedTuple{name_variables, NTuple{1, Vector{Float64}}}(tuple(eachcol(u0_NamedTuple[:SWATI][:,:,1])...)),
            GWAT   = NamedTuple{name_variables,      NTuple{1, Float64}}(u0_NamedTuple[:GWAT][:,:,1]),
            INTS   = NamedTuple{name_variables,      NTuple{1, Float64}}(u0_NamedTuple[:INTS][:,:,1]),
            INTR   = NamedTuple{name_variables,      NTuple{1, Float64}}(u0_NamedTuple[:INTR][:,:,1]),
            SNOW   = NamedTuple{name_variables,      NTuple{1, Float64}}(u0_NamedTuple[:SNOW][:,:,1]),
            RWU    = NamedTuple{name_variables,      NTuple{1, Float64}}(u0_NamedTuple[:RWU]        ),
            XYLEM  = NamedTuple{name_variables,      NTuple{1, Float64}}(u0_NamedTuple[:XYLEM]      ),
            CC     = NamedTuple{name_variables[[1]], NTuple{1, Float64}}(u0_NamedTuple[:CC][:,:,1]),
            SNOWLQ = NamedTuple{name_variables[[1]], NTuple{1, Float64}}(u0_NamedTuple[:SNOWLQ][:,:,1]),

            TRANI = NamedTuple{name_variables, NTuple{1, Vector{Float64}}}(tuple(eachcol(u0_NamedTuple[:TRANI][:,:,1])...)),

            aux   = NamedTuple{name_aux, NTuple{3, Vector{Float64}}}(tuple(eachcol(u0_NamedTuple[:aux])...)),
            accum = NamedTuple{name_accum, NTuple{31, Float64}}((0. for i in eachindex(name_accum))))
    end

    # # Give ArrayPartition as u0 to DiffEq.jl
    # u0 = ArrayPartition(u0_NamedTuple...)
    # u0_field_names = keys(u0_NamedTuple) # and save names of u0 to parameter vector

    # return
    return u0
end

function init_LWFB90_u0!(;u0::ComponentArray, continuous_SPAC, soil_discr, p_soil)

    N_iso = ifelse(continuous_SPAC.solver_options.simulate_isotopes, 2, 0)

    # A) Define initial conditions of states
    u_SWATIinit_mm      = LWFBrook90.KPT.FTheta(LWFBrook90.KPT.FWETNES(soil_discr["PSIM_init"], p_soil), p_soil) .*
                          p_soil.p_SWATMAX ./ p_soil.p_THSAT # see l.2020: https://github.com/pschmidtwalter/LWFBrook90R/blob/6f23dc1f6be9e1723b8df5b188804da5acc92e0f/src/md_brook90.f95#L2020

    u0.GWAT   .= continuous_SPAC.continuousIC.scalar[1:(N_iso+1), "u_GWAT_init_mm"]
    u0.INTS   .= continuous_SPAC.continuousIC.scalar[1:(N_iso+1), "u_INTS_init_mm"]
    u0.INTR   .= continuous_SPAC.continuousIC.scalar[1:(N_iso+1), "u_INTR_init_mm"]
    u0.SNOW   .= continuous_SPAC.continuousIC.scalar[1:(N_iso+1), "u_SNOW_init_mm"]
    u0.CC     .= continuous_SPAC.continuousIC.scalar[1,           "u_CC_init_MJ_per_m2"]
    u0.SNOWLQ .= continuous_SPAC.continuousIC.scalar[1,           "u_SNOWLQ_init_mm"]
    u0.SWATI.mm  .= u_SWATIinit_mm
    if (N_iso == 2)
        u0.SWATI.d18O  .= soil_discr["d18O_init"]
        u0.SWATI.d2H  .= soil_discr["d2H_init"]
    end

    u0.RWU.mm   = 0
    u0.XYLEM.mm = 5
    u0.TRANI.mm  = zeros(soil_discr["NLAYER"])
    if (N_iso == 2)
        u0.RWU.d18O   = soil_discr["d18O_init"][1] # start out with same concentration as in first soil layer
        u0.RWU.d2H    = soil_discr["d2H_init"][1]   # start out with same concentration as in first soil layer
        u0.XYLEM.d18O = soil_discr["d18O_init"][1] # start out with same concentration as in first soil layer
        u0.XYLEM.d2H  = soil_discr["d2H_init"][1]  # start out with same concentration as in first soil layer
        u0.TRANI.d18O .= soil_discr["d18O_init"]    # start out with same concentration as in       soil layer
        u0.TRANI.d2H  .= soil_discr["d2H_init"]     # start out with same concentration as in       soil layer
    end
    # # TODO(bernhard): if species-specific uptakes add here a totalRWU *PER SPECIES*
    # # TODO(bernhard): if species-specific uptakes add here a Xylem value *PER SPECIES*
    # # TODO(bernhard): if species-specific uptakes add here an uptake vector *PER SPECIES*

    # B) Define initial conditions of auxiliary soil states
    u0.aux # TODO: θ, ψ, K

    # C) Define initial conditions of accumulation variables
    if continuous_SPAC.solver_options.compute_intermediate_quantities
        # initialize terms for balance errors with initial values from u0
        new_SWAT       = sum(u_SWATIinit_mm) # total soil water in all layers, mm
        new_totalWATER = u0.INTR.mm + u0.INTS.mm + u0.SNOW.mm + new_SWAT + u0.GWAT.mm # total water in all compartments, mm

        u0.accum.totalSWAT = new_SWAT
        u0.accum.new_totalWATER = new_totalWATER
        # accums[:totalSWAT]       .= new_SWAT
        # accums[:new_totalWATER]  .= new_totalWATER
    end

    return nothing
end