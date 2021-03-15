module GLBLDECL

using Interpolations: interpolate, BSpline, Constant, Previous, scale, extrapolate, NoInterp
using DataFrames: DataFrame
using DataFramesMeta
using Dates: DateTime

using ..WAT: LWFRootGrowth

export derive_params_from_inputData

#######################
#######################
#######################
# Exported functions

"""
    derive_params_from_inputData(input_meteoveg::DataFrame, input_param::DataFrame, input_siteparam::DataFrame,
    input_precdat::DataFrame, input_pdur::DataFrame, input_soil_materials::DataFrame, input_soil_nodes::DataFrame)

Take input data sets and defines parameters needed for simulation.
"""
function derive_params_from_inputData(input_meteoveg::DataFrame,
                                  input_param::DataFrame,
                                  input_siteparam::DataFrame,
                                  input_precdat::DataFrame,
                                  input_pdur::DataFrame,
                                  input_soil_materials::DataFrame,
                                  input_soil_nodes::DataFrame,
                                  input_reference_date::DateTime)

    pfile_param = derive_params_from_input_param(input_param)
    # Defines: pfile_param[:ILAYER], pfile_param[:NDAYS], pfile_param[:HEAT], pfile_param[:ESLOPE], pfile_param[:ASPECT], pfile_param[:ALB], pfile_param[:ALBSN], pfile_param[:C1], pfile_param[:C2], pfile_param[:C3], pfile_param[:WNDRAT], pfile_param[:FETCH], pfile_param[:Z0W], pfile_param[:ZW], pfile_param[:LWIDTH], pfile_param[:Z0G], pfile_param[:Z0S], pfile_param[:LPC], pfile_param[:CS], pfile_param[:CZS], pfile_param[:CZR], pfile_param[:HS], pfile_param[:HR], pfile_param[:ZMINH], pfile_param[:RHOTP], pfile_param[:NN], pfile_param[:RSTEMP], pfile_param[:INTR_init], pfile_param[:INTS_init], pfile_param[:FRINTL], pfile_param[:FSINTL], pfile_param[:FRINTS], pfile_param[:FSINTS], pfile_param[:CINTRL], pfile_param[:CINTRS], pfile_param[:CINTSL], pfile_param[:CINTSS], pfile_param[:MELFAC], pfile_param[:CCFAC], pfile_param[:LAIMLT], pfile_param[:SAIMLT], pfile_param[:GRDMLT], pfile_param[:MAXLQF], pfile_param[:KSNVP], pfile_param[:SNODEN], pfile_param[:GLMAX], pfile_param[:CR], pfile_param[:GLMIN], pfile_param[:RM], pfile_param[:R5], pfile_param[:CVPD], pfile_param[:TL], pfile_param[:T1], pfile_param[:T2], pfile_param[:TH], pfile_param[:MXKPL], pfile_param[:MXRTLN], pfile_param[:inirlen], pfile_param[:inirdep], pfile_param[:rgroper], pfile_param[:FXYLEM], pfile_param[:PSICR], pfile_param[:RTRAD], pfile_param[:NOOUTF], pfile_param[:FXYLEM], pfile_param[:inirlen], pfile_param[:NLAYER], pfile_param[:ILAYER], pfile_param[:QLAYER], pfile_param[:IMODEL], pfile_param[:RSSA], pfile_param[:RSSB], pfile_param[:INFEXP], pfile_param[:BYPAR], pfile_param[:QFPAR], pfile_param[:QFFC], pfile_param[:IMPERV], pfile_param[:DSLOPE], pfile_param[:LENGTH], pfile_param[:DRAIN], pfile_param[:GSC], pfile_param[:GSP], pfile_param[:DTIMAX], pfile_param[:DSWMAX], pfile_param[:DPSIMX], )
    pfile_siteparam = derive_params_from_input_siteparam(input_siteparam)
    # Defines pfile_siteparam["p_LAT"]
    pfile_precdat = Nothing # derive_params_from_input_precdat(input_precdat)
    # Defines pfile_precdat
    pfile_pdur = derive_params_from_input_pdur(input_pdur)
    # Defines: pfile_pdur["DURATN"]
    pfile_soil = derive_params_from_input_soil(input_soil_materials,
                                               input_soil_nodes,
                                               pfile_param[:IMODEL],#IMODEL,
                                               pfile_param[:ILAYER],#ILAYER,
                                               pfile_param[:QLAYER],#QLAYER,
                                               pfile_param[:NLAYER],#NLAYER,
                                               pfile_param[:HEAT],#HEAT,
                                               pfile_param[:nmat],#nmat)
                                               pfile_param[:inirdep],
                                               pfile_param[:rgrorate])
    # Defines: pfile_soil["THICK"],pfile_soil["PSIM_init"],pfile_soil["frelden"],pfile_soil["PAR"],pfile_soil["STONEF"],pfile_soil["tini"],pfile_soil["HeatCapOld"],pfile_soil["TopInfT"],

    pfile_meteoveg = derive_params_from_input_meteoveg(input_meteoveg, input_reference_date,
                                                # for precipitation:
                                                pfile_siteparam["p_NPINT"],
                                                # for RootGrowth in LWFBrook90.jl:
                                                pfile_param[:NLAYER],
                                                pfile_param[:inirdep],
                                                pfile_param[:inirlen],
                                                pfile_param[:rgroper],
                                                pfile_soil["tini"],
                                                pfile_soil["frelden"])
    # Defines: pfile_meteoveg["p_GLOBRAD"], pfile_meteoveg["p_TMAX"], pfile_meteoveg["p_TMIN"],
    #          pfile_meteoveg["p_VAPPRES"], pfile_meteoveg["p_WIND"], pfile_meteoveg["p_PREC"], pfile_meteoveg["p_MESFL"],
    #          pfile_meteoveg["p_DENSEF"], pfile_meteoveg["p_HEIGHT"], pfile_meteoveg["p_LAI"], pfile_meteoveg["p_SAI"], pfile_meteoveg["p_AGE"],
    #          pfile_meteoveg["p_FRELDEN"]

    # TODO(bernhard): document input parameters: inirdep, inirlen, rgroper, tini, frelden,
    #                 as they are only used here, but not documented in comments in
    #                 func_DiffEq_definition_p.jl

    return (pfile_meteoveg, pfile_param, pfile_siteparam, pfile_precdat, pfile_pdur, pfile_soil)
end


#######################
#######################
#######################
# Unexported functions
"""
    derive_params_from_input_meteoveg(input_meteoveg::DataFrame)

Take climate and vegetation parameters in `input_meteoveg` and generates continuous parameters.
"""
function derive_params_from_input_meteoveg(
    input_meteoveg::DataFrame,
    input_reference_date::DateTime,
    # for precipitation
    p_NPINT,
    # root density parameters
    NLAYER,
    p_inirdep,
    p_inirlen,
    p_rgroper,
    p_tini,
    p_frelden)

    # 2) Interpolate input data in time
    time_range = range(minimum(input_meteoveg.days), maximum(input_meteoveg.days), length=length(input_meteoveg.days))

    # using Plots
    # time_range = range(minimum(input_meteoveg.days), maximum(input_meteoveg.days), length=length(input_meteoveg.days))
    # ts = 0:0.01:365
    # scatter(input_meteoveg.days, input_meteoveg.PRECIN)
    # plot!(ts, scale(interpolate(input_meteoveg.PRECIN, (BSpline(Constant{Previous}()))), time_range)(ts), label = "PRECIN {Previous}",  xlims=(0,30))
    # plot!(ts, scale(interpolate(input_meteoveg.PRECIN, (BSpline(Constant{Next}()))),     time_range)(ts), label = "PRECIN {Next}",      xlims=(0,30))
    # plot!(ts, scale(interpolate(input_meteoveg.PRECIN, (BSpline(Constant()))),           time_range)(ts), label = "PRECIN {Nearest}",   xlims=(0,30))
    # # plot!(ts, scale(interpolate(input_meteoveg.PRECIN, (BSpline(Linear()))),             time_range)(ts), label = "PRECIN Linear",   xlims=(0,30))

    p_GLOBRAD = extrapolate(scale(interpolate(input_meteoveg.GLOBRAD, (BSpline(Constant{Previous}()))), time_range) ,0)
    p_TMAX    = extrapolate(scale(interpolate(input_meteoveg.TMAX,    (BSpline(Constant{Previous}()))), time_range) ,0)
    p_TMIN    = extrapolate(scale(interpolate(input_meteoveg.TMIN,    (BSpline(Constant{Previous}()))), time_range) ,0)
    p_VAPPRES = extrapolate(scale(interpolate(input_meteoveg.VAPPRES, (BSpline(Constant{Previous}()))), time_range) ,0)
    p_WIND    = extrapolate(scale(interpolate(input_meteoveg.WIND,    (BSpline(Constant{Previous}()))), time_range) ,0)
    p_MESFL   = extrapolate(scale(interpolate(input_meteoveg.MESFL,   (BSpline(Constant{Previous}()))), time_range) ,0)
    p_DENSEF  = extrapolate(scale(interpolate(input_meteoveg.DENSEF,  (BSpline(Constant{Previous}()))), time_range) ,0)
    p_HEIGHT  = extrapolate(scale(interpolate(input_meteoveg.HEIGHT,  (BSpline(Constant{Previous}()))), time_range) ,0)
    p_LAI     = extrapolate(scale(interpolate(input_meteoveg.LAI,     (BSpline(Constant{Previous}()))), time_range) ,0)
    p_SAI     = extrapolate(scale(interpolate(input_meteoveg.SAI,     (BSpline(Constant{Previous}()))), time_range) ,0)
    p_AGE     = extrapolate(scale(interpolate(input_meteoveg.AGE,     (BSpline(Constant{Previous}()))), time_range) ,0)
    p_PREC    = extrapolate(scale(interpolate(input_meteoveg.PRECIN,  (BSpline(Constant{Previous}()))), time_range) ,0)

    # 2a Compute time dependent root density parameters
    # Which is a vector quantity that is dependent on time:
    p_RELDEN_2Darray = fill(NaN, nrow(input_meteoveg), NLAYER)
    for i in 1:nrow(input_meteoveg)
        p_RELDEN_2Darray[i,:] = LWFRootGrowth(p_frelden, p_tini, input_meteoveg.AGE[i], p_rgroper, p_inirdep, p_inirlen, NLAYER)
    end
    p_RELDEN =  extrapolate(scale(interpolate(p_RELDEN_2Darray, (BSpline(Constant{Previous}()), NoInterp()) ),# 1st dimension: ..., 2nd dimension NoInterp()
                            time_range, 1:size(p_RELDEN_2Darray,2)),
                    0) # extrapolate with fillvalue = 0

    # 2b Compute precipitation rate for precipitation intervals
    # NOTE: PRECIN is already in mm/day from the input data set
    #       No transformation is needed for p_PREC.
    #       Just define p_DTP.
    if p_NPINT == 1
        p_DTP = 1 / p_NPINT
    else
        error("Case with multiple precipitation intervals (using PRECDAT) is not implemented.")
    end

    # #### TODO(bernhard): get rid of this part below
    # if (isequal(p_NPINT, 1))
    #     # One precipitation interval per day. I.e. length of precipitation interval:
    #     p_DTP = 1 # (days), time step for precipitation interval (for PRECIN and MESFL)
    #               # in this case it is equal to simulation interval p_DT = 1.0 day

    #     # BROOK90 computed:
    #     p_fT_PREINT = p_PREC(integrator.t) / p_DTP # (mm per interval), precipitation
    #     # FB: however, this seems not correct in cases where DTP ≠ 1. I'd rather compute:
    #     p_fT_PREINT = p_PREC(integrator.t) * p_DTP # (mm per interval), precipitation
    # else
    #     # multiple precipitation intervals per day: (not implemented)
    #     error("Case with multiple precipitation intervals (using PRECDAT) is not implemented.")
    #     p_DTP = 1/p_NPINT # (days), time step for precipitation interval (for PRECIN and MESFL)
    #     # This would result in use of e.g. hourly PRECIN and hourly MESFL

    #     # PRECIN (mm/day)
    #     # PREINT (mm per interval) = PRECIN*p_DTP
    #     error("Case where input file PRECDAT is used is not implemented.
    #            Reading PRECDAT should result in PREINT (precipitation amount per interval)")
    # end
    # #### END TODO

    return Dict([("p_GLOBRAD",p_GLOBRAD),
                 ("p_TMAX",p_TMAX),
                 ("p_TMIN",p_TMIN),
                 ("p_VAPPRES",p_VAPPRES),
                 ("p_WIND",p_WIND),
                 ("p_PREC",p_PREC),
                 ("p_DTP",p_DTP),
                 ("p_NPINT",p_NPINT),
                 ("p_MESFL",p_MESFL),
                 ("p_DENSEF",p_DENSEF),
                 ("p_HEIGHT",p_HEIGHT),
                 ("p_LAI",p_LAI),
                 ("p_SAI",p_SAI),
                 ("p_AGE",p_AGE),
                 ("p_RELDEN",p_RELDEN),
                 ("input_reference_date", input_reference_date)])
end

"""
    derive_params_from_input_pdur(input_pdur)

Define constant parameters from input_pdur. TO BE REDEFINED
"""
function derive_params_from_input_pdur(input_pdur)
    # Interception parameters -------
    DURATN = input_pdur[:,1]
    return Dict([("DURATN",DURATN)])
end

"""
    derive_params_from_input_siteparam(input_siteparam)

Define constant parameters from input_siteparam. TO BE REDEFINED
"""
function derive_params_from_input_siteparam(input_siteparam)
    p_LAT  = input_siteparam[1,3]/57.296

    u_GWAT_init = input_siteparam[1,5]
    u_SNOW_init = input_siteparam[1,4]
    p_NPINT     = input_siteparam[1,6] # was precip_interval

    return Dict([("p_LAT",p_LAT),
                 ("u_GWAT_init",u_GWAT_init),
                 ("u_SNOW_init",u_SNOW_init),
                 ("p_NPINT",p_NPINT)])
end

"""
    derive_params_from_input_soil(input_soil_materials,IMODEL, ILAYER, QLAYER, NLAYER)

Define constant parameters from input_pdur. TO BE REDEFINED
"""
function derive_params_from_input_soil(input_soil_materials, input_soil_nodes, IMODEL, ILAYER, QLAYER, NLAYER, HEAT, nmat, inirdep, rgrorate)
    # Parse soil parameter according to material at different depth given by
    # arguments:
    # input_soil_materials: A matrix of the 8 soil materials parameters.
    #                 When imodel = 1 (Mualem-van Genuchten), these refer to:
    #                       mat, ths, thr, alpha (m-1), npar, ksat (mm d-1), tort (-), stonef (-).
    #                 When imodel = 0 (Clapp-Hornberger):
    #                       mat, thsat, thetaf, psif (kPa), bexp, kf (mm d-1), wetinf (-), stonef (-).


    #if (NLAYER > ML) || (ILAYER > NLAYER) || (QLAYER > NLAYER)
    if (ILAYER > NLAYER) || (QLAYER > NLAYER)
        error("Failure of QLAYER and ILAYER < NLAYER < ML")
    end

    # Parse the parameters for each material depending on wheter we use iModel ==1 or ==2
    ParMat    = fill(NaN, (nmat, 10))
    StonefMat = fill(NaN, (nmat, 1))
    for i = 1:nmat
        if IMODEL == 0
        # Clapp-Hornberger
        #             input_soil_materials[i,1]  # mat
        ParMat[i,1] = input_soil_materials[i,2]  # θs
        ParMat[i,2] = input_soil_materials[i,3]  # θf
        ParMat[i,4] = input_soil_materials[i,4]  # ψf (kPa)
        ParMat[i,9] = input_soil_materials[i,5]  # bexp (-)
        ParMat[i,3] = input_soil_materials[i,6]  # kf (mm d-1)
        ParMat[i,10] = input_soil_materials[i,7] # wetinf (-)
        StonefMat[i,1] = input_soil_materials[i,8] # stonef (-)
        ParMat[i,5] = 0.
        ParMat[i,6] = 0.
        ParMat[i,7] = 0.
        ParMat[i,8] = 0.
        # ParMat[i,CH]> (θs, θf, kf, ψf, 0, 0, 0, 0, bexp, wetinf)
        end

        if IMODEL == 1
            # Mualem-van Genuchten
            #             input_soil_materials[i,1]  # mat
            ParMat[i,1] = input_soil_materials[i,2]  # θs
            ParMat[i,10] = input_soil_materials[i,3] # θr
            ParMat[i,7] = input_soil_materials[i,4]  # α (m-1)
            ParMat[i,8] = input_soil_materials[i,5]  # npar
            ParMat[i,6] = input_soil_materials[i,6]  # ksat (mm d-1)
            ParMat[i,9] = input_soil_materials[i,7]  # tort (-)
            StonefMat[i,1] = input_soil_materials[i,8] # stonef (-)
            ParMat[i,2] = 0.
            ParMat[i,4] = 0.
            ParMat[i,5] = 0.
            # ..it's a sin......
            # Hard default [mm d-1] for saturated hydraulic conductivity at field capacity
            ParMat[i,3] = 2            # K(θ_fc) = 2 (mm d-1)
            # ParMat_MvG ==> (θs,  0, K(θ_fc), 0, 0, ksat, α, npar, tort, θr)
        end
        # ParMat_CH  ==> (θs, θf, kf,     ψf, 0, 0,    0,    0, bexp, wetinf)
        # ParMat_MvG ==> (θs,  0, K(θ_fc), 0, 0, ksat, α, npar, tort, θr)
    end
    dep       = fill(NaN, NLAYER) # soil depth [m]  #TODO(bernhard): not exported
    THICK     = fill(NaN, NLAYER)
    mat       = fill(0, NLAYER)   # material_id for each soil layer #TODO(bernhard): not exported
    PSIM_init = fill(NaN, NLAYER)
    frelden   = fill(NaN, NLAYER)
    for i = 1:NLAYER
        if HEAT == 1
            dep[i]       = input_soil_nodes[i,2]
            THICK[i]     = input_soil_nodes[i,3]
            mat[i]       = Int( input_soil_nodes[i,4] )
            PSIM_init[i] = input_soil_nodes[i,5]
            frelden[i]   = input_soil_nodes[i,6]
            # TemperatureNew(i)       = input_soil_nodes[i,7] we don't have it in the input file!!!
            if i > NLAYER
                MUE[i] = THICK[i] / ( THICK[i] + THICK(I+1) )
                ZL[i]  = 0.5 * ( THICK[i] + THICK(I+1) )
            else
                MUE[i] = 0.5
                ZL[i]  = THICK[i]
            end
            TMean[i]   = 0.
        else
            dep[i]       = input_soil_nodes[i,2]
            THICK[i]     = input_soil_nodes[i,3]
            mat[i]       = Int( input_soil_nodes[i,4] )
            PSIM_init[i] = input_soil_nodes[i,5]
            frelden[i]   = input_soil_nodes[i,6]
        end
    end
    depmax = dep[1] - THICK[1] / 1000.

    #..from material-specific to layer-specific parameter values
    PAR    = fill(NaN, (NLAYER, 10)) # MPAR=10
    STONEF = fill(NaN, (NLAYER))
    for i = 1:NLAYER
        for j = 1:10 # MPAR=10
            PAR[i,j] = ParMat[mat[i], j]
            # PAR_CH  ==> (θs, θf, kf,     ψf, 0, 0,    0,    0, bexp, wetinf)
            # PAR_MvG ==> (θs,  0, K(θ_fc), 0, 0, ksat, α, npar, tort, θr)
        end
        STONEF[i] = StonefMat[ mat[i] ]
    end
    if IMODEL == 0
        # Clapp-Hornberger
        PAR = DataFrame(PAR, ["θs","θf","kf","ψf","dummy1","dummy2","dummy3","dummy4","bexp","wetinf"])
    elseif IMODEL == 1
        # Mualem-van Genuchten
        PAR = DataFrame(PAR, ["θs","dummy1","K(θ_fc)","dummy2","dummy3","Ksat","α","n","tort","θr"])
    end

    # find thickness of maximum root zone
    # frelden: relative values of final root density per unit volume
    i1 = findfirst(frelden .> 1.e-6) # first layer where frelden[i] is >=1.e-6
    i2 = findlast(frelden .> 1.e-6)  # last layer (in 1:NLAYER) where frelden[i] is >=1.e-6
    if !isnothing(i1) && isnothing(i2)
        i2 = NLAYER
    end

    tini = fill(NaN, NLAYER)
    for i = 1:NLAYER
        tini[i] = 1.e+20 # initial time for root growth in layer
        if i >= i1 && i <= i2
            frelden[i] = max( frelden[i], 1.01e-6)
        end
        if frelden[i] >= 1.e-6 && (depmax-dep[i]) <= inirdep
            tini[i] = 0.
        end
        if frelden[i] >= 1.e-6 && (depmax-dep[i]) > inirdep
            if rgrorate > 0
                tini[i] = (depmax-dep[i]-inirdep)/rgrorate
            end
        end
        # write(*,*)'dep= ',dep[i],' tini= ',tini[i]
    end


    # heat flow -------
    # we assign some so compilation does not complain
    TopInfT = 1
    #       if (HEAT == 1)
    #        READ (12,*) Comment
    #        READ (12,*) tTop, Comment
    #        tTop=tTop
    #        READ (12,*) tBot, Comment
    #        tBot=tBot
    #        READ (12,*) TopInfT, Comment
    #        READ (12,*) BotInfT, Comment
    #        READ (12,*) kTopT, Comment
    #        READ (12,*) kBotT, Comment
    #        DO 207 I = 1, 7
    #         READ (12,*) Comment
    # 207    CONTINUE
    #        DO 208 I = 1, nmat
    #         READ (12,*) ilay, SV[i], OV[i], HB1[i], HB2[i], HB3[i]
    #         TPar[1,I) = SV[i]
    #         TPar[2,I) = OV[i]
    #         TPar[3,I) = THDis
    # C        thermal conductivities -- transfer from [J m-1 s-1 K-1] to  [J mm-1 d-1 K-1]
    #         TPar[4,I) = HB1[i] * 86.400
    #         TPar[5,I) = HB2[i] * 86.400
    #         TPar[6,I) = HB3[i] * 86.400
    # C         volumetric heat capacities for solid, water and organic -- transfer from [MJ m-2 mm-1 K-1] to [J mm-3 K-1]
    #         TPar[7,I) = p_CVSOL   # To define in module_CONSTANTS.jl: p_CVSOL = ?  # CVSOL  - volumetric heat capacity of solid (MJ m-2 mm-1 K-1)
    #         TPar[8,I) = p_CVORG   # To define in module_CONSTANTS.jl: p_CVORG = ?  # volumetric heat capacity of organic material (MJ m-2 mm-1 K-1) (hillel98)
    #         TPar[9,I) = LWFBrook90.CONSTANTS.p_CVLQ
    # 208    CONTINUE
    #        READ (12,*) C
    #       end

    ### # from LWFBrook90R:md_brook90.f95 (lines 105)
    TPar = fill(NaN, (nmat, 10))
    TPar .= -1
    HeatCapOld = fill(NaN, NLAYER)
    for i = 1:NLAYER
        HeatCapOld[i] = TPar[mat[i],7] * TPar[mat[i],1] + TPar[mat[i],8]
                    # TPar[2,mat[i])+TPar[9,mat[i])*SWATI[i]/THICK[i]
    end
    ###



    return Dict([
                ("THICK",THICK),
                ("PSIM_init",PSIM_init),
                ("frelden",frelden),
                ("PAR",PAR),
                ("STONEF",STONEF),
                ("tini",tini),
                #
                ("HeatCapOld",HeatCapOld),
                ("TopInfT", TopInfT)])
end

"""
    derive_params_from_input_param(input_param)

Defines constant parameters from input_param. TO BE REDEFINED
"""
function derive_params_from_input_param(input_param)
    return copy(input_param[1,:])
end


end