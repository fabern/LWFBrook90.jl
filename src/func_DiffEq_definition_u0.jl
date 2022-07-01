"""
    define_LWFB90_u0()

Generate vector u0 needed for ODE() problem in DiffEq.jl package.
"""
function define_LWFB90_u0(p, uScalar_initial,
    uSoil_initial_ψM_kPa, uSoil_initial_δ18O_permil, uSoil_initial_δ2H_permil,
    compute_intermediate_quantities;
    simulate_isotopes::Bool = false)

    numericNaN = -9999.99 # instead of NaN in the u-vector set them to this value

    # A) Define initial conditions

    # A1) amounts:
    u_GWAT_init_mm      = uScalar_initial[1, "u_GWAT_init_mm"]
    u_INTS_init_mm      = uScalar_initial[1, "u_INTS_init_mm"]
    u_INTR_init_mm      = uScalar_initial[1, "u_INTR_init_mm"]
    u_SNOW_init_mm      = uScalar_initial[1, "u_SNOW_init_mm"]
    u_CC_init_MJ_per_m2 = uScalar_initial[1, "u_CC_init_MJ_per_m2"]
    u_SNOWLQ_init_mm    = uScalar_initial[1, "u_SNOWLQ_init_mm"]
    p_soil = p[1][1]
    # u_SWATIinit_mm    = LWFBrook90.KPT.FWETNES(uSoil_initial_ψM_kPa, p_soil) .*
    #                      p_soil.p_SWATMAX
    u_SWATIinit_mm      = LWFBrook90.KPT.FTheta(LWFBrook90.KPT.FWETNES(uSoil_initial_ψM_kPa, p_soil), p_soil) .*
                          p_soil.p_SWATMAX ./ p_soil.p_THSAT # see l.2020: https://github.com/pschmidtwalter/LWFBrook90R/blob/6f23dc1f6be9e1723b8df5b188804da5acc92e0f/src/md_brook90.f95#L2020

    N_separate_treespecies = 1 # TODO(bernhard) currently only one species is implemented
    u_totalRWUinit_mmday= [0 for j in 1:N_separate_treespecies]
    u_Xyleminit_mm      = [5 for j in 1:N_separate_treespecies] # mm, Volume of well mixed xylem storage per area, TODO(bernhard): make this a parameter
    u_TRANIinit_mmday   = vcat([zeros(size(u_SWATIinit_mm)) for j in 1:N_separate_treespecies]...)
    # TODO(bernhard): if species-specific uptakes add here a totalRWU *PER SPECIES*
    # TODO(bernhard): if species-specific uptakes add here a Xylem value *PER SPECIES*
    # TODO(bernhard): if species-specific uptakes add here an uptake vector *PER SPECIES*

    # A2) isotopic concentrations:
    if simulate_isotopes
        # isotopes d18O
        u_GWAT_init_d18O      = uScalar_initial[2, "u_GWAT_init_mm"]
        u_INTS_init_d18O      = uScalar_initial[2, "u_INTS_init_mm"]
        u_INTR_init_d18O      = uScalar_initial[2, "u_INTR_init_mm"]
        u_SNOW_init_d18O      = uScalar_initial[2, "u_SNOW_init_mm"]
        # u_CC_init_MJ_per_m2_d18O = uScalar_initial[2, "u_CC_init_MJ_per_m2"] # state has no isotopic signature
        # u_SNOWLQ_init_d18O       = uScalar_initial[2, "u_SNOWLQ_init_mm"]    # state has no isotopic signature
        u_SWATIinit_δ18O      = uSoil_initial_δ18O_permil
        u_RWUinit_d18O        = [uSoil_initial_δ18O_permil[1] for j in 1:N_separate_treespecies]   # start out with same concentration as in first soil layer
        u_XYLEMinit_d18O      = [uSoil_initial_δ18O_permil[1] for j in 1:N_separate_treespecies]   # start out with same concentration as in first soil layer
        u_TRANIinit_d18O      = vcat([u_SWATIinit_δ18O for j in 1:N_separate_treespecies]...)      # start out with same concentration as all soil layers

        # isotopes d2H
        u_GWAT_init_d2H      = uScalar_initial[3, "u_GWAT_init_mm"]
        u_INTS_init_d2H      = uScalar_initial[3, "u_INTS_init_mm"]
        u_INTR_init_d2H      = uScalar_initial[3, "u_INTR_init_mm"]
        u_SNOW_init_d2H      = uScalar_initial[3, "u_SNOW_init_mm"]
        # u_CC_init_MJ_per_m2_d2H = uScalar_initial[3, "u_CC_init_MJ_per_m2"] # state has no isotopic signature
        # u_SNOWLQ_init_d2H       = uScalar_initial[3, "u_SNOWLQ_init_mm"]    # state has no isotopic signature
        u_SWATIinit_δ2H      = uSoil_initial_δ2H_permil
        u_RWUinit_d2H        = [uSoil_initial_δ2H_permil[1] for j in 1:N_separate_treespecies]   # start out with same concentration as in first soil layer
        u_XYLEMinit_d2H      = [uSoil_initial_δ2H_permil[1] for j in 1:N_separate_treespecies]   # start out with same concentration as in first soil layer
        u_TRANIinit_d2H      = vcat([u_SWATIinit_δ2H for j in 1:N_separate_treespecies]...)      # start out with same concentration as all soil layers
    end

    # A3) accumulation variables:
    # initialize terms for balance errors with initial values from u0
    new_SWAT       = sum(u_SWATIinit_mm) # total soil water in all layers, mm
    new_totalWATER = u_INTR_init_mm + u_INTS_init_mm + u_SNOW_init_mm + new_SWAT + u_GWAT_init_mm # total water in all compartments, mm

    # B) build state vector (actually state array)
    # B1) simple case:
    u0 = [u_GWAT_init_mm          ;
          u_INTS_init_mm          ;
          u_INTR_init_mm          ;
          u_SNOW_init_mm          ;
          u_CC_init_MJ_per_m2     ;
          u_SNOWLQ_init_mm        ;
          u_SWATIinit_mm...       ;
          u_totalRWUinit_mmday... ; # TODO(bernhard): if species-specific uptakes add here a totalRWU *PER SPECIES*
          u_Xyleminit_mm...       ; # TODO(bernhard): if species-specific uptakes add here a Xylem value *PER SPECIES*
          u_TRANIinit_mmday...    ] # TODO(bernhard): if species-specific uptakes add here an uptake vector *PER SPECIES*

    # B2) isotopic concentrations:
    if simulate_isotopes
        u0 = [u_GWAT_init_mm           u_GWAT_init_d18O           u_GWAT_init_d2H;
              u_INTS_init_mm           u_INTS_init_d18O           u_INTS_init_d2H;
              u_INTR_init_mm           u_INTR_init_d18O           u_INTR_init_d2H;
              u_SNOW_init_mm           u_SNOW_init_d18O           u_SNOW_init_d2H;
              u_CC_init_MJ_per_m2      numericNaN                 numericNaN;
              u_SNOWLQ_init_mm         numericNaN                 numericNaN;
              hcat(u_SWATIinit_mm,     u_SWATIinit_δ18O,          u_SWATIinit_δ2H);
              hcat(u_totalRWUinit_mmday, u_RWUinit_d18O,          u_RWUinit_d2H);    # TODO(bernhard): if species-specific uptakes add here a totalRWU *PER SPECIES*
              hcat(u_Xyleminit_mm,       u_XYLEMinit_d18O,        u_XYLEMinit_d2H);  # TODO(bernhard): if species-specific uptakes add here a Xylem value *PER SPECIES*
              hcat(u_TRANIinit_mmday,    u_TRANIinit_d18O,        u_TRANIinit_d2H);  # TODO(bernhard): if species-specific uptakes add here an uptake vector *PER SPECIES*
        ]
    end

    # B3) accumulation variables:
    if compute_intermediate_quantities
        N_accum_var = 31
        u_accum_init = fill(0., N_accum_var)
        u_accum_init[28] = new_SWAT
        u_accum_init[29] = new_totalWATER

        if simulate_isotopes
            u0 = [u_GWAT_init_mm        u_GWAT_init_d18O        u_GWAT_init_d2H;
                  u_INTS_init_mm        u_INTS_init_d18O        u_INTS_init_d2H;
                  u_INTR_init_mm        u_INTR_init_d18O        u_INTR_init_d2H;
                  u_SNOW_init_mm        u_SNOW_init_d18O        u_SNOW_init_d2H;
                  u_CC_init_MJ_per_m2   numericNaN              numericNaN;
                  u_SNOWLQ_init_mm      numericNaN              numericNaN;
                  hcat(u_SWATIinit_mm,        u_SWATIinit_δ18O,             u_SWATIinit_δ2H);
                  hcat(u_totalRWUinit_mmday,  u_RWUinit_d18O,               u_RWUinit_d2H);   # TODO(bernhard): if species-specific uptakes add here a totalRWU *PER SPECIES*
                  hcat(u_Xyleminit_mm,        u_XYLEMinit_d18O,             u_XYLEMinit_d2H); # TODO(bernhard): if species-specific uptakes add here a Xylem value *PER SPECIES*
                  hcat(u_TRANIinit_mmday,     u_TRANIinit_d18O,             u_TRANIinit_d2H); # TODO(bernhard): if species-specific uptakes add here an uptake vector *PER SPECIES*
                  hcat(u_accum_init,          fill(numericNaN,N_accum_var), fill(numericNaN,N_accum_var));
            ]
        else
            u0 = [u_GWAT_init_mm          ;
                  u_INTS_init_mm          ;
                  u_INTR_init_mm          ;
                  u_SNOW_init_mm          ;
                  u_CC_init_MJ_per_m2     ;
                  u_SNOWLQ_init_mm        ;
                  u_SWATIinit_mm...       ;
                  u_totalRWUinit_mmday... ; # TODO(bernhard): if species-specific uptakes add here a totalRWU *PER SPECIES*
                  u_Xyleminit_mm...       ; # TODO(bernhard): if species-specific uptakes add here a Xylem value *PER SPECIES*
                  u_TRANIinit_mmday...    ; # TODO(bernhard): if species-specific uptakes add here an uptake vector *PER SPECIES*
                  u_accum_init]
        end
    end

    # C) Define how to access the values:
    row_idx_scalars = (GWAT     = 1 ,
                       INTS     = 2 ,
                       INTR     = 3 ,
                       SNOW     = 4 ,
                       CC       = 5 ,
                       SNOWLQ   = 6 ,
                       totalRWU = 48,
                       XylemV   = 49,)
    row_idx_SWATI = collect(6 + layer_idx                 for layer_idx in 1:p_soil.NLAYER)
    row_idx_RWU   = collect(8 + layer_idx + p_soil.NLAYER for layer_idx in 1:p_soil.NLAYER)
    if compute_intermediate_quantities
        row_idx_accum = length(row_idx_scalars) + length(row_idx_SWATI) +length(row_idx_RWU) .+ (1:N_accum_var)
        names_accum = reshape([
                "cum_d_prec";"cum_d_rfal";"cum_d_sfal";"cum_d_rint"; "cum_d_sint";"cum_d_rsno";
                "cum_d_rnet";"cum_d_smlt";"cum_d_evap";"cum_d_tran";"cum_d_irvp";"cum_d_isvp";
                "cum_d_slvp";"cum_d_snvp";"cum_d_pint";"cum_d_ptran";"cum_d_pslvp";
                "flow";"seep";"srfl";"slfl";"byfl";"dsfl";"gwfl";"vrfln";
                "cum_d_rthr";"cum_d_sthr";
                "totalSWAT"; "new_totalWATER"; "BALERD_SWAT"; "BALERD_total";
            ],1,:)
    else
        row_idx_accum = ([])
        names_accum = []
    end
    # idx_u_scalar_amounts       -> row_idx_scalars, 1
    # idx_u_scalar_isotopes_d18O -> row_idx_scalars, 2
    # idx_u_scalar_isotopes_d2H  -> row_idx_scalars, 3
    # idx_u_vector_amounts       -> row_idx_SWATI, 1
    # idx_u_vector_isotopes_d18O -> row_idx_SWATI, 2
    # idx_u_vector_isotopes_d2H  -> row_idx_SWATI, 3


    # write the access patterns into the parameter vector p
    # (note, because a tuple is immutable we need to generate a new tuple p)
    p_cst = (p[1][1],
        p[1][2],
        p[1][3],
            (FLAG_MualVanGen = p[1][4][:FLAG_MualVanGen],
            compute_intermediate_quantities = p[1][4][:compute_intermediate_quantities],
            simulate_isotopes = p[1][4][:simulate_isotopes],

            row_idx_scalars = row_idx_scalars,
            row_idx_SWATI   = row_idx_SWATI,
            row_idx_RWU     = row_idx_RWU,
            row_idx_accum   = row_idx_accum,
            names_accum     = names_accum,
            col_idx_d18O    = 2,
            col_idx_d2H     = 3))
    p_fT = p[2]
    p_fu = p[3]
    p_cache = p[4]
    p = (p_cst, p_fT, p_fu, p_cache)

    return u0, p
end