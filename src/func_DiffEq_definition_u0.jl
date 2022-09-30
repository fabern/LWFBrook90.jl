"""
    define_LWFB90_u0()

Generate vector u0 needed for ODE() problem in DiffEq.jl package.
"""
# simulate_isotopes               = continuous_SPAC.solver_options.simulate_isotopes
# compute_intermediate_quantities = continuous_SPAC.solver_options.compute_intermediate_quantities
# NLAYER = 21
# u0, u0_field_names, u0_variable_names, names_accum = define_LWFB90_u0(;simulate_isotopes, compute_intermediate_quantities, NLAYER)
# p = define_LWFB90_p(continuous_SPAC, soil_discr)
# p_soil = p[1][1]
# init_LWFB90_u0!(;u0=u0, u0_field_names=u0_field_names, u0_variable_names=u0_variable_names,
#                  continuous_SPAC=continuous_SPAC, soil_discr=soil_discr, p_soil=p[1][1])
function define_LWFB90_u0(;simulate_isotopes, compute_intermediate_quantities, NLAYER)
    u0_variable_names = simulate_isotopes ? (d18O = 2, d2H = 3) : ()
    # u0_treespecies_names = (default = 1) # e.g. (beech = 1, oak = 2, spruce = 3)
    names_accum = compute_intermediate_quantities ? reshape([
                "cum_d_prec";"cum_d_rfal";"cum_d_sfal";"cum_d_rint"; "cum_d_sint";"cum_d_rsno";
                "cum_d_rnet";"cum_d_smlt";"cum_d_evap";"cum_d_tran";"cum_d_irvp";"cum_d_isvp";
                "cum_d_slvp";"cum_d_snvp";"cum_d_pint";"cum_d_ptran";"cum_d_pslvp";
                "flow";"seep";"srfl";"slfl";"byfl";"dsfl";"gwfl";"vrfln";
                "cum_d_rthr";"cum_d_sthr";
                "totalSWAT"; "new_totalWATER"; "BALERD_SWAT"; "BALERD_total";
            ],1,:) : []

    N_isotopes             = length(u0_variable_names)
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
    # Give ArrayPartition as u0 to DiffEq.jl
    u0 = ArrayPartition(u0_NamedTuple...)
    u0_field_names = keys(u0_NamedTuple) # and save names of u0 to parameter vector
    # In f() or cb(): get back a NamedTuple for easy access/modification
    # states = NamedTuple{p_u0Names}(u0_AP.x)
    # # which allows then to modify it as
    # u0_AP.x[1]
    # states.GWAT[1] = 0.33
    # states.GWAT[2] = -13
    # states.GWAT[3] = -95
    # u0_AP.x[1]

    # return
    return u0, u0_field_names, u0_variable_names, names_accum
end

function init_LWFB90_u0!(;u0::ArrayPartition, u0_field_names, u0_variable_names,
                         continuous_SPAC, soil_discr, p_soil) #names_accum, simulate_isotopes, compute_intermediate_quantities)

    # simulate_isotopes               = continuous_SPAC.solver_options.simulate_isotopes
    # compute_intermediate_quantities = continuous_SPAC.solver_options.compute_intermediate_quantities

    states = NamedTuple{u0_field_names}(u0.x)
    # accums = NamedTuple{Tuple(Symbol.(names_accum[:]))}(states.accum)

    N_iso = length(u0_variable_names)

    # A) Define initial conditions of states
    u_SWATIinit_mm      = LWFBrook90.KPT.FTheta(LWFBrook90.KPT.FWETNES(soil_discr["PSIM_init"], p_soil), p_soil) .*
                          p_soil.p_SWATMAX ./ p_soil.p_THSAT # see l.2020: https://github.com/pschmidtwalter/LWFBrook90R/blob/6f23dc1f6be9e1723b8df5b188804da5acc92e0f/src/md_brook90.f95#L2020

    states.GWAT   .= continuous_SPAC.continuousIC.scalar[1:(N_iso+1), "u_GWAT_init_mm"]'
    states.INTS   .= continuous_SPAC.continuousIC.scalar[1:(N_iso+1), "u_INTS_init_mm"]'
    states.INTR   .= continuous_SPAC.continuousIC.scalar[1:(N_iso+1), "u_INTR_init_mm"]'
    states.SNOW   .= continuous_SPAC.continuousIC.scalar[1:(N_iso+1), "u_SNOW_init_mm"]'
    states.CC     .= continuous_SPAC.continuousIC.scalar[1,           "u_CC_init_MJ_per_m2"]'
    states.SNOWLQ .= continuous_SPAC.continuousIC.scalar[1,           "u_SNOWLQ_init_mm"]'
    states.SWATI  .= [u_SWATIinit_mm soil_discr["d18O_init"] soil_discr["d2H_init"]]

    states.RWU   .= [0                           soil_discr["d18O_init"][1] soil_discr["d2H_init"][1]] # start out with same concentration as in first soil layer
    states.XYLEM .= [5                           soil_discr["d18O_init"][1] soil_discr["d2H_init"][1]] # start out with same concentration as in first soil layer
    states.TRANI .= [zeros(soil_discr["NLAYER"]) soil_discr["d18O_init"]    soil_discr["d2H_init"]   ] # start out with same concentration as in       soil layer
    # # TODO(bernhard): if species-specific uptakes add here a totalRWU *PER SPECIES*
    # # TODO(bernhard): if species-specific uptakes add here a Xylem value *PER SPECIES*
    # # TODO(bernhard): if species-specific uptakes add here an uptake vector *PER SPECIES*
    # TODO: as we're not using u0_variable_names, we are implicitly assuming the order: mm, d18O, d2H
    @assert u0_variable_names == (d18O = 2, d2H = 3)

    # B) Define initial conditions of auxiliary soil states
    states.aux # TODO: θ, ψ, K

    # C) Define initial conditions of accumulation variables
    # initialize terms for balance errors with initial values from u0
    new_SWAT       = sum(u_SWATIinit_mm) # total soil water in all layers, mm
    new_totalWATER = states.INTR[1] + states.INTS[1] + states.SNOW[1] + new_SWAT + states.GWAT[1] # total water in all compartments, mm

    states.accum[28] = new_SWAT
    states.accum[29] = new_totalWATER
    # accums[:totalSWAT]       .= new_SWAT
    # accums[:new_totalWATER]  .= new_totalWATER

    return nothing
end