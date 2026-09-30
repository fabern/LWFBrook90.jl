"""
    define_LWFB90_u0()

Generate vector u0 needed for ODE() problem in DiffEq.jl package.
"""
function define_LWFB90_u0(;simulate_isotopes, compute_intermediate_quantities, NLAYER, species_names::Vector{Symbol} = [:species1])
    name_states = ifelse(simulate_isotopes, (:mm,    :d18O, :d2H), (:mm,))
    name_fluxes = ifelse(simulate_isotopes, (:mmday, :d18O, :d2H), (:mmday,))
    name_aux       = (:θ,:ψ,:K)
    name_accum     = (:cum_d_prec, :cum_d_rfal, :cum_d_sfal, :cum_d_rint,  :cum_d_sint, :cum_d_rsno,
                    :cum_d_rnet, :cum_d_smlt, :evap, :cum_d_tran, :cum_d_irvp, :cum_d_isvp,
                    :cum_d_slvp, :cum_d_snvp, :cum_d_pint, :cum_d_ptran, :cum_d_pslvp,
                    :flow, :seep, :srfl, :slfl, :byfl, :dsfl, :gwfl, :vrfln,
                    :cum_d_rthr, :cum_d_sthr, :cum_d_irrig,
                    :StorageSWAT,  :StorageWATER,  :BALERD_SWAT,  :BALERD_total)

    if !isempty(species_names) && species_names != [:species1]
        name_accum = (name_accum..., (Symbol("$(name)_$sp") for sp in species_names
                                      for name in (:cum_d_tran, :cum_d_irvp, :cum_d_isvp, :cum_d_ptran, :cum_d_pint))...)
    end

    variable_names = simulate_isotopes ? (d18O = 2, d2H = 3) : ()
    N_isotopes             = length(variable_names)
    N_separate_treespecies = 1 # Legacy aggregate states use one column; species states are added below # TODO: check renaming to clarify this. This refers to N_aggregated (which is by obviously 1)
    N_accum_var            = length(name_accum)

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

    # Add one set of interception and hydraulic states per species.
    species_components = Pair{Symbol, Any}[]
    if !isempty(species_names) && species_names != [:species1]
        for sp in species_names
            push!(species_components, Symbol("INTS_$sp") => NamedTuple{name_states, NTuple{1+N_isotopes, Float64}}(tuple(zeros(1+N_isotopes)...)))
            push!(species_components, Symbol("INTR_$sp") => NamedTuple{name_states, NTuple{1+N_isotopes, Float64}}(tuple(zeros(1+N_isotopes)...)))
            push!(species_components, Symbol("RWU_$sp") => NamedTuple{name_fluxes, NTuple{1+N_isotopes, Float64}}(tuple(zeros(1+N_isotopes)...)))
            push!(species_components, Symbol("XYLEM_$sp") => NamedTuple{name_states, NTuple{1+N_isotopes, Float64}}(tuple(zeros(1+N_isotopes)...)))
            push!(species_components, Symbol("TRANI_$sp") => NamedTuple{name_fluxes, NTuple{1+N_isotopes, Vector{Float64}}}(tuple(fill(zeros(NLAYER), 1+N_isotopes)...)))
        end
    end

    # Give ComponentArray as u0 to DiffEq.jl
    if simulate_isotopes
        u0 = ComponentArray(
            SWATI  = NamedTuple{name_states, NTuple{3, Vector{Float64}}}(tuple(eachcol(u0_NamedTuple[:SWATI][:,:,1])...)),
            GWAT   = NamedTuple{name_states,      NTuple{3, Float64}}(u0_NamedTuple[:GWAT][:,:,1]),
            INTS   = NamedTuple{name_states,      NTuple{3, Float64}}(u0_NamedTuple[:INTS][:,:,1]),
            INTR   = NamedTuple{name_states,      NTuple{3, Float64}}(u0_NamedTuple[:INTR][:,:,1]),
            SNOW   = NamedTuple{name_states,      NTuple{3, Float64}}(u0_NamedTuple[:SNOW][:,:,1]),
            RWU    = NamedTuple{name_fluxes,      NTuple{3, Float64}}(u0_NamedTuple[:RWU]        ),
            XYLEM  = NamedTuple{name_states,      NTuple{3, Float64}}(u0_NamedTuple[:XYLEM]      ),
            CC     = NamedTuple{(:MJm2,),         NTuple{1, Float64}}(u0_NamedTuple[:CC][:,:,1]),
            SNOWLQ = NamedTuple{name_states[[1]], NTuple{1, Float64}}(u0_NamedTuple[:SNOWLQ][:,:,1]),

            TRANI = NamedTuple{name_fluxes, NTuple{3, Vector{Float64}}}(tuple(eachcol(u0_NamedTuple[:TRANI][:,:,1])...)),

            aux   = NamedTuple{name_aux,   NTuple{3, Vector{Float64}}}(tuple(eachcol(u0_NamedTuple[:aux])...)),
            accum = NamedTuple{name_accum, NTuple{N_accum_var, Float64}}((0. for i in eachindex(name_accum)));
            species_components...)
    else
        # TODO(bernhard): check if this is bad programming if NTuple{1, ...} depends on runtime variable simulate_isotopes...
        u0 = ComponentArray(
            SWATI  = NamedTuple{name_states, NTuple{1, Vector{Float64}}}(tuple(eachcol(u0_NamedTuple[:SWATI][:,:,1])...)),
            GWAT   = NamedTuple{name_states,      NTuple{1, Float64}}(u0_NamedTuple[:GWAT][:,:,1]),
            INTS   = NamedTuple{name_states,      NTuple{1, Float64}}(u0_NamedTuple[:INTS][:,:,1]),
            INTR   = NamedTuple{name_states,      NTuple{1, Float64}}(u0_NamedTuple[:INTR][:,:,1]),
            SNOW   = NamedTuple{name_states,      NTuple{1, Float64}}(u0_NamedTuple[:SNOW][:,:,1]),
            RWU    = NamedTuple{name_fluxes,      NTuple{1, Float64}}(u0_NamedTuple[:RWU]        ),
            XYLEM  = NamedTuple{name_states,      NTuple{1, Float64}}(u0_NamedTuple[:XYLEM]      ),
            CC     = NamedTuple{(:MJm2,),         NTuple{1, Float64}}(u0_NamedTuple[:CC][:,:,1]),
            SNOWLQ = NamedTuple{name_states[[1]], NTuple{1, Float64}}(u0_NamedTuple[:SNOWLQ][:,:,1]),

            TRANI = NamedTuple{name_fluxes, NTuple{1, Vector{Float64}}}(tuple(eachcol(u0_NamedTuple[:TRANI][:,:,1])...)),

            aux   = NamedTuple{name_aux, NTuple{3, Vector{Float64}}}(tuple(eachcol(u0_NamedTuple[:aux])...)),
            accum = NamedTuple{name_accum, NTuple{N_accum_var, Float64}}((0. for i in eachindex(name_accum)));
            species_components...)
    end

    # # Give ArrayPartition as u0 to DiffEq.jl
    # u0 = ArrayPartition(u0_NamedTuple...)
    # u0_field_names = keys(u0_NamedTuple) # and save names of u0 to parameter vector

    # return
    return u0
end

function init_LWFB90_u0!(;u0::ComponentArray, parametrizedSPAC, p_soil)

    N_iso = ifelse(parametrizedSPAC.solver_options.simulate_isotopes, 2, 0)
    species_names = parametrizedSPAC.pars.species_names
    cover_fractions = parametrizedSPAC.pars.cover_fractions

    soil_PSIM_init = parametrizedSPAC.soil_discretization.df.uAux_PSIM_init_kPa
    soil_d18O_init = parametrizedSPAC.soil_discretization.df.u_delta18O_init_permil
    soil_d2H_init  = parametrizedSPAC.soil_discretization.df.u_delta2H_init_permil

    # A) Define initial conditions of states
    u_SWATIinit_mm      = LWFBrook90.KPT.FTheta(LWFBrook90.KPT.FWETNES(soil_PSIM_init, p_soil), p_soil) .*
                          p_soil.p_SWATMAX ./ p_soil.p_THSAT # see l.2020: https://github.com/pschmidtwalter/LWFBrook90R/blob/6f23dc1f6be9e1723b8df5b188804da5acc92e0f/src/md_brook90.f95#L2020

    u0.GWAT   .= parametrizedSPAC.pars.IC_scalar[1:(N_iso+1), "u_GWAT_init_mm"]
    u0.INTS   .= parametrizedSPAC.pars.IC_scalar[1:(N_iso+1), "u_INTS_init_mm"]
    u0.INTR   .= parametrizedSPAC.pars.IC_scalar[1:(N_iso+1), "u_INTR_init_mm"]
    u0.SNOW   .= parametrizedSPAC.pars.IC_scalar[1:(N_iso+1), "u_SNOW_init_mm"]
    u0.CC     .= parametrizedSPAC.pars.IC_scalar[1,           "u_CC_init_MJ_per_m2"]
    u0.SNOWLQ .= parametrizedSPAC.pars.IC_scalar[1,           "u_SNOWLQ_init_mm"]
    u0.SWATI.mm  .= u_SWATIinit_mm
    if (N_iso == 2)
        u0.SWATI.d18O  .= soil_d18O_init
        u0.SWATI.d2H  .= soil_d2H_init
    end

    if !isempty(species_names) && species_names != [:species1]
        sum_ints = 0.0
        sum_intr = 0.0
        for sp in species_names
            w_s = cover_fractions[sp]
            ints_col = "u_INTS_init_mm_$sp"
            intr_col = "u_INTR_init_mm_$sp"

            sp_ints = hasproperty(parametrizedSPAC.pars.IC_scalar, Symbol(ints_col)) ? parametrizedSPAC.pars.IC_scalar[1:(N_iso+1), ints_col] : (parametrizedSPAC.pars.IC_scalar[1:(N_iso+1), "u_INTS_init_mm"] .* [w_s; ones(N_iso)])
            sp_intr = hasproperty(parametrizedSPAC.pars.IC_scalar, Symbol(intr_col)) ? parametrizedSPAC.pars.IC_scalar[1:(N_iso+1), intr_col] : (parametrizedSPAC.pars.IC_scalar[1:(N_iso+1), "u_INTR_init_mm"] .* [w_s; ones(N_iso)])

            u0[Symbol("INTS_$sp")] .= sp_ints
            u0[Symbol("INTR_$sp")] .= sp_intr
            sum_ints += sp_ints[1]
            sum_intr += sp_intr[1]

            getproperty(u0, Symbol("RWU_$sp")).mmday = 0.0
            getproperty(u0, Symbol("XYLEM_$sp")).mm = 5.0
            getproperty(u0, Symbol("TRANI_$sp")).mmday .= zeros(nrow(parametrizedSPAC.soil_discretization.df))

        end
        u0.INTS.mm = sum_ints
        u0.INTR.mm = sum_intr
        if N_iso == 2
            u0.INTS.d18O = u0[Symbol("INTS_$(species_names[1])")].d18O
            u0.INTS.d2H  = u0[Symbol("INTS_$(species_names[1])")].d2H
            u0.INTR.d18O = u0[Symbol("INTR_$(species_names[1])")].d18O
            u0.INTR.d2H  = u0[Symbol("INTR_$(species_names[1])")].d2H
        end
    end

    u0.RWU.mmday   = 0
    u0.XYLEM.mm    = 5
    u0.TRANI.mmday = zeros(nrow(parametrizedSPAC.soil_discretization.df))
    if (N_iso == 2)
        u0.RWU.d18O   = soil_d18O_init[1] # start out with same concentration as in first soil layer
        u0.RWU.d2H    = soil_d2H_init[1]   # start out with same concentration as in first soil layer
        u0.XYLEM.d18O = soil_d18O_init[1] # start out with same concentration as in first soil layer
        u0.XYLEM.d2H  = soil_d2H_init[1]  # start out with same concentration as in first soil layer
        u0.TRANI.d18O .= soil_d18O_init    # start out with same concentration as in       soil layer
        u0.TRANI.d2H  .= soil_d2H_init     # start out with same concentration as in       soil layer
    end

    # Species-specific uptake, xylem, and layer uptake states are initialized above.

    # B) Define initial conditions of auxiliary soil states
    u0.aux # TODO: θ, ψ, K

    # C) Define initial conditions of accumulation variables
    if parametrizedSPAC.solver_options.compute_intermediate_quantities
        # initialize terms for balance errors with initial values from u0
        new_SWAT       = sum(u_SWATIinit_mm) # total soil water in all layers, mm
        StorageWATER = u0.INTR.mm + u0.INTS.mm + u0.SNOW.mm + new_SWAT + u0.GWAT.mm # total water in all compartments, mm

        u0.accum.StorageSWAT = new_SWAT
        u0.accum.StorageWATER = StorageWATER
        # accums[:StorageSWAT]       .= new_SWAT
        # accums[:StorageWATER]  .= StorageWATER
    end

    return nothing
end