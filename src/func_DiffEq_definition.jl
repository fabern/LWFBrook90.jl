""" define_diff_eq_parameters(NLAYER, IMODEL, constant_dt_solver, NOOUTF, Reset,
    pfile_meteo, pfile_siteparam, pfile_param, pfile_soil, pfile_pdur)\n 
    Generates vector p needed for ODE() probelm in DiffEq.jl package."""
function define_DiffEq_parameters(NLAYER, IMODEL, constant_dt_solver, NOOUTF, Reset,
    pfile_meteo, pfile_siteparam, pfile_param, pfile_soil, pfile_pdur)

    # 

    ########
    # 1) Parse pfile inputs:

    # heat flow (unimplemented)
    p_HEAT    = pfile_param["HEAT"]
    p_TopInfT = pfile_soil["TopInfT"]

    # radiation aspect parameters
    p_LAT   = pfile_siteparam["p_LAT"]
    p_GLMAX = pfile_param["GLMAX"]
    p_GLMIN = pfile_param["GLMIN"]
    p_ESLOPE= pfile_param["ESLOPE"]
    p_DSLOPE= pfile_param["DSLOPE"]
    # equivalent slope for radiation calculations
    p_L1, p_L2 = LWFBrook90Julia.SUN.EQUIVSLP(p_LAT, p_ESLOPE, pfile_param["ASPECT"])

    p_SNODEN = pfile_param["SNODEN"]
    p_MXRTLN = pfile_param["MXRTLN"]
    p_MXKPL  = pfile_param["MXKPL"]
    p_CS     = pfile_param["CS"]
    p_Z0S    = pfile_param["Z0S"]
    p_Z0G    = pfile_param["Z0G"]
    p_ZMINH  = pfile_param["ZMINH"]
    p_CZS    = pfile_param["CZS"]
    p_CZR    = pfile_param["CZR"]
    p_HS     = pfile_param["HS"]
    p_HR     = pfile_param["HR"]
    p_LPC    = pfile_param["LPC"]
    p_RTRAD  = pfile_param["RTRAD"]
    p_FXYLEM = pfile_param["FXYLEM"]
    p_WNDRAT = pfile_param["WNDRAT"]
    p_FETCH  = pfile_param["FETCH"]
    p_Z0W    = pfile_param["Z0W"]
    p_ZW     = pfile_param["ZW"]
    p_RSTEMP = pfile_param["RSTEMP"]
    p_LWIDTH = pfile_param["LWIDTH"]
    p_RHOTP  = pfile_param["RHOTP"]
    p_NN     = pfile_param["NN"]
    p_KSNVP  = pfile_param["KSNVP"]
    p_ALBSN  = pfile_param["ALBSN"]
    p_ALB    = pfile_param["ALB"]
    p_RSSA   = pfile_param["RSSA"]
    p_RSSB   = pfile_param["RSSB"]
    p_CCFAC  = pfile_param["CCFAC"]
    p_MELFAC = pfile_param["MELFAC"]
    p_LAIMLT = pfile_param["LAIMLT"]
    p_SAIMLT = pfile_param["SAIMLT"]
    p_GRDMLT = pfile_param["GRDMLT"]
    p_C1     = pfile_param["C1"]
    p_C2     = pfile_param["C2"]
    p_C3     = pfile_param["C3"]
    p_CR     = pfile_param["CR"]
    p_R5     = pfile_param["R5"]
    p_CVPD   = pfile_param["CVPD"]
    p_RM     = pfile_param["RM"]
    p_TL     = pfile_param["TL"]
    p_T1     = pfile_param["T1"]
    p_T2     = pfile_param["T2"]
    p_TH     = pfile_param["TH"]
    p_FSINTL = pfile_param["FSINTL"]
    p_FSINTS = pfile_param["FSINTS"]
    p_CINTSL = pfile_param["CINTSL"]
    p_CINTSS = pfile_param["CINTSS"]
    p_FRINTL = pfile_param["FRINTL"]
    p_FRINTS = pfile_param["FRINTS"]
    p_CINTRL = pfile_param["CINTRL"]
    p_CINTRS = pfile_param["CINTRS"]
    p_MAXLQF = pfile_param["MAXLQF"]
    p_QFPAR  = pfile_param["QFPAR"]
    p_QFFC   = pfile_param["QFFC"]
    p_IMPERV = pfile_param["IMPERV"]
    p_BYPAR  = pfile_param["BYPAR"]
    p_LENGTH = pfile_param["LENGTH"]
    p_DPSIMX = pfile_param["DPSIMX"]
    p_DRAIN  = pfile_param["DRAIN"]
    p_DTIMAX = pfile_param["DTIMAX"]
    p_DSWMAX = pfile_param["DSWMAX"]
    p_GSC    = pfile_param["GSC"]
    p_GSP    = pfile_param["GSP"]

    p_DURATN = pfile_pdur["DURATN"]

    # unused: # root density parameters
    # unused: p_inirdep= pfile_param["inirdep"] # for RootGrowth in LWFBrook90Julia
    # unused: p_inirlen= pfile_param["inirlen"] # for RootGrowth in LWFBrook90Julia
    # unused: p_rgroper= pfile_param["rgroper"] # for RootGrowth in LWFBrook90Julia
    # unused: p_tini   = pfile_soil["tini"]     # for RootGrowth in LWFBrook90Julia
    # unused: p_frelden= pfile_soil["frelden"]  # for RootGrowth in LWFBrook90Julia
    # unused: p_HeatCapOld = pfile_soil["HeatCapOld"]

    # soil water parameters
    # unused u_aux_PSIMinit = pfile_soil["PSIM_init"]
    p_PSICR= pfile_param["PSICR"]
    p_THICK= pfile_soil["THICK"]
    p_STONEF=pfile_soil["STONEF"]
    p_THSAT =pfile_soil["PAR"][!,"θs"]
    if IMODEL == 0
        p_THETAF = pfile_soil["PAR"][!,"θf"]
        p_KF     = pfile_soil["PAR"][!,"kf"]
        p_PSIF   = pfile_soil["PAR"][!,"ψf"]
        p_BEXP   = pfile_soil["PAR"][!,"bexp"]
        p_WETINF = pfile_soil["PAR"][!,"wetinf"]
        p_Kθfc   = fill(NaN, NLAYER)
        p_Ksat   = fill(NaN, NLAYER)
        p_MvGα   = fill(NaN, NLAYER)
        p_MvGn   = fill(NaN, NLAYER)
        p_MvGl   = fill(NaN, NLAYER)
        p_θr     = fill(NaN, NLAYER)
    elseif IMODEL == 1
        p_THETAF = fill(NaN, NLAYER)
        p_KF     = fill(NaN, NLAYER)
        p_PSIF   = fill(NaN, NLAYER)
        p_BEXP   = fill(NaN, NLAYER)
        p_WETINF = fill(NaN, NLAYER)
        p_Kθfc = pfile_soil["PAR"][!,"K(θ_fc)"]
        p_Ksat = pfile_soil["PAR"][!,"Ksat"]
        p_MvGα = pfile_soil["PAR"][!,"α"]
        p_MvGn = pfile_soil["PAR"][!,"n"]
        p_MvGl = pfile_soil["PAR"][!,"tort"]
        p_θr   = pfile_soil["PAR"][!,"θr"]
    else
        error("Unsupported IMODEL: $IMODEL")
    end

    # infiltration parameters INFPAR()
    p_ILAYER = pfile_param["ILAYER"] # number of layers over which infiltration is distributed
    p_INFEXP = pfile_param["INFEXP"]
    p_INFRAC = LWFBrook90Julia.WAT.INFPAR(p_INFEXP, p_ILAYER, p_THICK, NLAYER)

    # soil water parameters
    (p_PSIG,             # gravity potential (kPa),
    p_SWATMX,           # maximum water storage for layer (mm),
    p_WETF,             # wetness at field capacity (dimensionless),
    #p_WETc, # wetness at field capacity (dimensionless) # TODO(bernhard): was unused with Brook90R, check if used in LWFBrook90R
    p_CHM,              # Clapp and Hornberger m (kPa),
    p_CHN,              # Clapp and Hornberger n,
    p_KSAT,             # saturated hydraulic conductivity (mm/d)
    p_PSIF,
    p_THETAF
    # TODO(bernhard): clean up use of p_PSIF and p_THETAF (depending on IMODEL they are 
    #                 initialized as NaN. Are they overwritten by SOILPAR?)
    ) = LWFBrook90Julia.KPT.SOILPAR(
            LWFBrook90Julia.CONSTANTS.p_RHOWG,
            p_THICK, # layer thicknesses (mm)
            p_THETAF,# θ at field capacity (-)
            p_THSAT, # θ at saturation == matrix porosity (-)

            p_STONEF,# stone volume fraction, unitless
            p_BEXP,  # exponent for psi-theta relation
            p_KF,    # hydraulic conductivity at field capacity (mm/d)
            p_PSIF,  # ψ at field capacity (kPa)
            p_WETINF,# wetness at dry end of near-saturation range

            p_Kθfc,
            p_PSICR, # minimum plant leaf water potential (MPa)
            p_Ksat, p_MvGl, p_MvGn, p_MvGα, p_θr,

            NLAYER, IMODEL)

    # source area parameters SRFPAR()
    p_QLAYER = # number of soil layers for SRFL
        pfile_param["QLAYER"]
    p_SWATQX = # maximum water storage for layers 1 through QLAYER (mm)
        sum(p_SWATMX[1:pfile_param["QLAYER"]])
    p_SWATQF = # water storage at field capacity for layers 1 through QLAYER (mm)
        sum(p_THETAF[1:pfile_param["QLAYER"]] .* p_THICK[1:pfile_param["QLAYER"]] .* (1 .- p_STONEF[1:pfile_param["QLAYER"]]))


    ########
    # 2) Define parameters for differential equation:
    # 2a) Constant parameters

    # p_cst_1 for both RHS and CallBack in DiffEq.jl
    p_cst_1 = (constant_dt_solver, NLAYER,
        p_SWATMX, p_PSIF, p_BEXP, p_WETINF, p_WETF, p_CHM, p_CHN, p_PSIG, p_KF,

        # FOR MSBITERATE:
        p_QLAYER, p_SWATQX, p_QFPAR, p_SWATQF, p_QFFC, p_IMPERV,
        p_LENGTH, p_DSLOPE, LWFBrook90Julia.CONSTANTS.p_RHOWG, p_DPSIMX,
        p_KSAT, p_DRAIN, p_DTIMAX, p_INFRAC, p_DSWMAX,
        p_GSC, p_GSP,
        
        # FOR UNIMPLEMENTED HEAT FLOW:
        p_HEAT, p_TopInfT,

        p_BYPAR)

    # p_cst_2 only for CallBack in DiffEq.jl
    p_cst_2 = (p_LAT, p_ESLOPE, p_L1, p_L2,
        p_SNODEN, p_MXRTLN, p_MXKPL, p_CS,
        p_Z0S, p_Z0G,
        p_ZMINH, p_CZS, p_CZR, p_HS, p_HR, p_LPC,
        p_THICK, p_STONEF, p_RTRAD, p_FXYLEM,
        p_WNDRAT, p_FETCH, p_Z0W, p_ZW,
        p_RSTEMP,
        LWFBrook90Julia.CONSTANTS.p_CVICE,
        p_LWIDTH, p_RHOTP, p_NN, p_KSNVP,
        p_ALBSN, p_ALB,
        p_RSSA, p_RSSB,
        p_CCFAC, p_MELFAC, p_LAIMLT, p_SAIMLT,

        LWFBrook90Julia.CONSTANTS.p_WTOMJ, p_C1, p_C2, p_C3, p_CR,
        p_GLMIN, p_GLMAX, p_R5, p_CVPD, p_RM, p_TL, p_T1, p_T2, p_TH,
        p_PSICR, NOOUTF,

        # for MSBPREINT:
        p_FSINTL, p_FSINTS, p_CINTSL, p_CINTSS,
        p_FRINTL, p_FRINTS, p_CINTRL, p_CINTRS,
        p_DURATN, p_MAXLQF, p_GRDMLT)

    p_cst = (p_cst_1, p_cst_2)
    
    # 2b) Time varying parameters (e.g. meteorological forcings)
    dt_input_data = 1 # 24*3600
    p_fT = (LWFBrook90Julia.p_DOY, LWFBrook90Julia.p_MONTHN, 
            pfile_meteo["p_GLOBRAD"], 
            pfile_meteo["p_TMAX"], 
            pfile_meteo["p_TMIN"], 
            pfile_meteo["p_VAPPRES"], 
            pfile_meteo["p_WIND"], 
            pfile_meteo["p_PRECIN"], 
            dt_input_data, 
            pfile_meteo["p_MESFL"],
            pfile_meteo["p_DENSEF"], 
            pfile_meteo["p_HEIGHT"], 
            pfile_meteo["p_LAI"], 
            pfile_meteo["p_SAI"], 
            pfile_meteo["p_AGE"])

    # dt_input_data = 1/24 # 1*3600
    # p_fT = (p_DOY, p_MONTHN, p_GLOBRAD, p_TMAX, p_TMIN, p_VAPPRES, p_WIND, p_PRECIN_h, dt_input_data, p_MESFL_h)

    # 2c) Time varying "parameters" (depending on state variables) 
    #     These need to be exchanged between CallBack and RHS in DiffEq.jl which is why they
    #     can temporarily be saved in the parameter vector to avoid computing them twice
    
    # Initialize placeholder for parameters that depend on solution and are computed
    p_fu = [NaN, NaN, fill(NaN, NLAYER), NaN] # Use array instead of tuple to be able to mutate

    # 3) Return different types of parameters as a single object
    return (p_cst, p_fT, p_fu)
end