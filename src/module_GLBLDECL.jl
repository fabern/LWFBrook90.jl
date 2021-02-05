# fabian.bernhard@wsl.ch, 2021-01-02

module GLBLDECL

using DataFrames
using Interpolations: interpolate, BSpline, Constant, scale, extrapolate
using DataFramesMeta

export define_Globals_ParametersAndVariables
export derive_params_from_input

#######################
#######################
#######################
# Exported functions

"""define_Globals_ParametersAndVariables()\n Function that defines many constant parameters. TO BE REDEFINED.
"""
function define_Globals_ParametersAndVariables()
    return Nothing
    # TODO(bernhard): is this needed?
end

"""
derive_params_from_input(input_meteo::DataFrame, input_param::DataFrame, input_siteparam::DataFrame,
input_precdat::DataFrame, input_pdur::DataFrame, input_soil_materials::DataFrame, input_soil_nodes::DataFrame)\n
Takes input data sets and defines parameters needed for simulation.
"""
function derive_params_from_input(input_meteo::DataFrame,
                                  input_param::DataFrame,
                                  input_siteparam::DataFrame,
                                  input_precdat::DataFrame,
                                  input_pdur::DataFrame,
                                  input_soil_materials::DataFrame,
                                  input_soil_nodes::DataFrame)

    pfile_meteo = derive_params_from_input_meteo(input_meteo)
    # Defines: pfile_meteo["p_GLOBRAD"], pfile_meteo["p_TMAX"], pfile_meteo["p_TMIN"], pfile_meteo["p_VAPPRES"], pfile_meteo["p_WIND"], pfile_meteo["p_PRECIN"], pfile_meteo["p_MESFL"], pfile_meteo["p_DENSEF"], pfile_meteo["p_HEIGHT"], pfile_meteo["p_LAI"], pfile_meteo["p_SAI"], pfile_meteo["p_AGE"]
    pfile_param = derive_params_from_input_param(input_param)
    # Defines: pfile_param["ILAYER"], pfile_param["NDAYS"], pfile_param["HEAT"], pfile_param["ESLOPE"], pfile_param["ASPECT"], pfile_param["ALB"], pfile_param["ALBSN"], pfile_param["C1"], pfile_param["C2"], pfile_param["C3"], pfile_param["WNDRAT"], pfile_param["FETCH"], pfile_param["Z0W"], pfile_param["ZW"], pfile_param["LWIDTH"], pfile_param["Z0G"], pfile_param["Z0S"], pfile_param["LPC"], pfile_param["CS"], pfile_param["CZS"], pfile_param["CZR"], pfile_param["HS"], pfile_param["HR"], pfile_param["ZMINH"], pfile_param["RHOTP"], pfile_param["NN"], pfile_param["RSTEMP"], pfile_param["INTR_init"], pfile_param["INTS_init"], pfile_param["FRINTL"], pfile_param["FSINTL"], pfile_param["FRINTS"], pfile_param["FSINTS"], pfile_param["CINTRL"], pfile_param["CINTRS"], pfile_param["CINTSL"], pfile_param["CINTSS"], pfile_param["MELFAC"], pfile_param["CCFAC"], pfile_param["LAIMLT"], pfile_param["SAIMLT"], pfile_param["GRDMLT"], pfile_param["MAXLQF"], pfile_param["KSNVP"], pfile_param["SNODEN"], pfile_param["GLMAX"], pfile_param["CR"], pfile_param["GLMIN"], pfile_param["RM"], pfile_param["R5"], pfile_param["CVPD"], pfile_param["TL"], pfile_param["T1"], pfile_param["T2"], pfile_param["TH"], pfile_param["MXKPL"], pfile_param["MXRTLN"], pfile_param["inirlen"], pfile_param["inirdep"], pfile_param["rgroper"], pfile_param["FXYLEM"], pfile_param["PSICR"], pfile_param["RTRAD"], pfile_param["NOOUTF"], pfile_param["FXYLEM"], pfile_param["inirlen"], pfile_param["NLAYER"], pfile_param["ILAYER"], pfile_param["QLAYER"], pfile_param["IMODEL"], pfile_param["RSSA"], pfile_param["RSSB"], pfile_param["INFEXP"], pfile_param["BYPAR"], pfile_param["QFPAR"], pfile_param["QFFC"], pfile_param["IMPERV"], pfile_param["DSLOPE"], pfile_param["LENGTH"], pfile_param["DRAIN"], pfile_param["GSC"], pfile_param["GSP"], pfile_param["DTIMAX"], pfile_param["DSWMAX"], pfile_param["DPSIMX"], )
    pfile_siteparam = derive_params_from_input_siteparam(input_siteparam)
    # Defines pfile_siteparam["p_LAT"]
    pfile_precdat = Nothing # derive_params_from_input_precdat(input_precdat)
    # Defines pfile_precdat
    pfile_pdur = derive_params_from_input_pdur(input_pdur)
    # Defines: pfile_pdur["DURATN"]
    pfile_soil = derive_params_from_input_soil(input_soil_materials, 
                                               input_soil_nodes, 
                                               pfile_param["IMODEL"],#IMODEL, 
                                               pfile_param["ILAYER"],#ILAYER, 
                                               pfile_param["QLAYER"],#QLAYER, 
                                               pfile_param["NLAYER"],#NLAYER, 
                                               pfile_param["HEAT"],#HEAT,
                                               pfile_param["nmat"],#nmat)
                                               pfile_param["inirdep"],
                                               pfile_param["rgrorate"])
    # Defines: pfile_soil["THICK"],pfile_soil["PSIM_init"],pfile_soil["frelden"],pfile_soil["PAR"],pfile_soil["STONEF"],pfile_soil["tini"],pfile_soil["HeatCapOld"],pfile_soil["TopInfT"], 
    
    return (pfile_meteo, pfile_param, pfile_siteparam, pfile_precdat, pfile_pdur, pfile_soil)
end


#######################
#######################
#######################
# Unexported functions
"""
derive_params_from_input_meteo(input_meteo::DataFrame)\n
Takes climate and vegetation parameters in `input_meteo` and generates continuous parameters.
"""
function derive_params_from_input_meteo(input_meteo::DataFrame)

    function interpolate_uniform(data, minScale, maxScale)
    @linq data |>
        #TODO(bernhard) make this work with option previous...
        interpolate(BSpline(Constant())) |>
        #variant b: interpolate(BSpline(Constant{Previous}())) |>
        #variant c: interpolate(BSpline(Constant{Next}())) |>
        # variant1 (assuming discrete scale with uniform spacing): scale(minScale:maxScale) |>
        # variant2 (assuming uniform spacing):
        scale(range(minScale, maxScale, length=length(data))) |>
        extrapolate(0) #fillvalue = 0

        # http://juliamath.github.io/Interpolations.jl/latest/control/#Scaled-BSplines-1
    end

    # 2) Interpolate input data in time
    p_GLOBRAD = interpolate_uniform(input_meteo[:,:GLOBRAD],minimum(input_meteo[:,:days]), maximum(input_meteo[:,:days]))
    p_TMAX    = interpolate_uniform(input_meteo[:,:TMAX],   minimum(input_meteo[:,:days]), maximum(input_meteo[:,:days]))
    p_TMIN    = interpolate_uniform(input_meteo[:,:TMIN],   minimum(input_meteo[:,:days]), maximum(input_meteo[:,:days]))
    p_VAPPRES = interpolate_uniform(input_meteo[:,:VAPPRES],minimum(input_meteo[:,:days]), maximum(input_meteo[:,:days]))
    p_WIND    = interpolate_uniform(input_meteo[:,:WIND],   minimum(input_meteo[:,:days]), maximum(input_meteo[:,:days]))
    p_PRECIN  = interpolate_uniform(input_meteo[:,:PRECIN], minimum(input_meteo[:,:days]), maximum(input_meteo[:,:days]))
    p_MESFL   = interpolate_uniform(input_meteo[:,:MESFL],  minimum(input_meteo[:,:days]), maximum(input_meteo[:,:days]))
    p_DENSEF  = interpolate_uniform(input_meteo[:,:DENSEF], minimum(input_meteo[:,:days]), maximum(input_meteo[:,:days]))
    p_HEIGHT  = interpolate_uniform(input_meteo[:,:HEIGHT], minimum(input_meteo[:,:days]), maximum(input_meteo[:,:days]))
    p_LAI     = interpolate_uniform(input_meteo[:,:LAI],    minimum(input_meteo[:,:days]), maximum(input_meteo[:,:days]))
    p_SAI     = interpolate_uniform(input_meteo[:,:SAI],    minimum(input_meteo[:,:days]), maximum(input_meteo[:,:days]))
    p_AGE     = interpolate_uniform(input_meteo[:,:AGE],    minimum(input_meteo[:,:days]), maximum(input_meteo[:,:days]))

    return Dict([("p_GLOBRAD",p_GLOBRAD),
                 ("p_TMAX",p_TMAX),
                 ("p_TMIN",p_TMIN),
                 ("p_VAPPRES",p_VAPPRES),
                 ("p_WIND",p_WIND),
                 ("p_PRECIN",p_PRECIN),
                 ("p_MESFL",p_MESFL),
                 ("p_DENSEF",p_DENSEF),
                 ("p_HEIGHT",p_HEIGHT),
                 ("p_LAI",p_LAI),
                 ("p_SAI",p_SAI),
                 ("p_AGE",p_AGE)])
end

"""derive_params_from_input_pdur(input_pdur)\n Function that defines constant parameters from input_pdur. TO BE REDEFINED
"""
function derive_params_from_input_pdur(input_pdur)
    # Interception parameters -------
    DURATN = input_pdur[:,1]
    return Dict([("DURATN",DURATN)])
end

"""derive_params_from_input_siteparam(input_siteparam)\n Function that defines constant parameters from input_siteparam. TO BE REDEFINED
"""
function derive_params_from_input_siteparam(input_siteparam)
    p_LAT  = input_siteparam[1,3]/57.296
    return Dict([("p_LAT",p_LAT)])
end

"""derive_params_from_input_soil(input_soil_materials, input_soil_nodes, IMODEL, ILAYER, QLAYER, NLAYER)\n Function that defines constant parameters from input_pdur. TO BE REDEFINED
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
    #         TPar[9,I) = p_CVLQ    # To define in module_CONSTANTS.jl: p_ThDis = ?  # Longitudinal thermal dispersivity (m)
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

"""derive_params_from_input_param(input_param)\nFunction that defines constant parameters from input_param. TO BE REDEFINED
"""
function derive_params_from_input_param(input_param)
    # from LWFBrook90R:PFILE.h


    # output specifications -------
    NDAYS = Int( input_param[ 1,1] )
    HEAT  = Int( input_param[ 2,1] ) # flag to include (1) or to exclude (0) soil heat flow

    # Meteorologic parameters -------
    ESLOPE = input_param[ 3,1]
    ASPECT = input_param[ 4,1]

    # Convert to radians
    ESLOPE = ESLOPE / 57.296
    ASPECT = ASPECT / 57.296

    ALB = input_param[ 5,1]
    ALBSN = input_param[ 6,1]
    C1 = input_param[ 7,1]
    C2 = input_param[ 8,1]
    C3 = input_param[ 9,1]
    WNDRAT = input_param[10,1]
    FETCH = input_param[11,1]
    Z0W = input_param[12,1]
    ZW = input_param[13,1]

    # Canopy parameters -------
    LWIDTH = input_param[14,1]
    Z0G = input_param[15,1]
    Z0S = input_param[16,1]
    LPC = input_param[17,1]
    CS = input_param[18,1]
    CZS = input_param[19,1]
    CZR = input_param[20,1]
    HS = input_param[21,1]
    HR = input_param[22,1]
    ZMINH = input_param[23,1]
    RHOTP = input_param[24,1]
    NN = input_param[25,1]

    # Interception parameters -------
    # DURATN = input_pdur[:,1]

    RSTEMP = input_param[26,1]
    INTR_init = input_param[27,1] # initial state
    INTS_init = input_param[28,1] # initial state
    FRINTL = input_param[29,1]
    FSINTL = input_param[30,1]
    FRINTS = input_param[31,1]
    FSINTS = input_param[32,1]
    CINTRL = input_param[33,1]
    CINTRS = input_param[34,1]
    CINTSL = input_param[35,1]
    CINTSS = input_param[36,1]

    # Snow parameters -------
    MELFAC = input_param[37,1]
    CCFAC = input_param[38,1]
    LAIMLT = input_param[39,1]
    SAIMLT = input_param[40,1]
    GRDMLT = input_param[41,1]
    MAXLQF = input_param[42,1]
    KSNVP = input_param[43,1]
    SNODEN = input_param[44,1]

    # leaf parameters affecting PE -------
    GLMAX = input_param[45,1]
    CR = input_param[46,1]
    GLMIN = input_param[47,1]
    RM = input_param[48,1]
    R5 = input_param[49,1]
    CVPD = input_param[50,1]
    TL = input_param[51,1]
    T1 = input_param[52,1]
    T2 = input_param[53,1]
    TH = input_param[54,1]

    # plant parameters affecting soil-water supply -------
    MXKPL = input_param[55,1]
    MXRTLN = input_param[56,1]
    inirlen = input_param[57,1]
    inirdep = input_param[58,1]  # only exported for derive_params_from_input_soil()
    rgrorate = input_param[59,1] # vertical root grow rate [m a-1] #TODO(bernhard): only exported for derive_params_from_input_soil()
    rgroper = input_param[60,1]
    FXYLEM = input_param[61,1]
    PSICR = input_param[62,1]
    RTRAD = input_param[63,1]
    NOOUTF = Int( input_param[64,1] )

    FXYLEM  = max(FXYLEM, 0.990)
    inirlen = min(inirlen, 0.010)
    inirdep = min(inirdep, 0.010)

    # soil parameters -------
    NLAYER = Int( input_param[65,1] )
    nmat = Int( input_param[66,1] )    # number of materials #TODO(bernhard): only exported for derive_params_from_input_soil()
    ILAYER = Int( input_param[67,1] )
    QLAYER = Int( input_param[68,1] )
    IMODEL = Int( input_param[69,1] )
    
    RSSA = input_param[70,1]
    RSSB = input_param[71,1]


    # flow parameters -------
    INFEXP = input_param[72,1]
    BYPAR = Int( input_param[73,1] )
    QFPAR = input_param[74,1]
    QFFC = input_param[75,1]
    IMPERV = input_param[76,1]
    DSLOPE = input_param[77,1]
    LENGTH = input_param[78,1]
    DRAIN = input_param[79,1]
    GSC = input_param[80,1]
    GSP = input_param[81,1]

    # integration parameters -------
    DTIMAX = input_param[82,1]
    DSWMAX = input_param[83,1]
    DPSIMX = input_param[84,1]

    return Dict([("NDAYS",NDAYS),
                ("HEAT",HEAT),
                ("ESLOPE",ESLOPE),
                ("ASPECT",ASPECT),
                ("ALB",ALB),
                ("ALBSN",ALBSN),
                ("C1",C1),
                ("C2",C2),
                ("C3",C3),
                ("WNDRAT",WNDRAT),
                ("FETCH",FETCH),
                ("Z0W",Z0W),
                ("ZW",ZW),
                ("LWIDTH",LWIDTH),
                ("Z0G",Z0G),
                ("Z0S",Z0S),
                ("LPC",LPC),
                ("CS",CS),
                ("CZS",CZS),
                ("CZR",CZR),
                ("HS",HS),
                ("HR",HR),
                ("ZMINH",ZMINH),
                ("RHOTP",RHOTP),
                ("NN",NN),
                ("RSTEMP",RSTEMP),
                ("INTR_init",INTR_init),
                ("INTS_init",INTS_init),
                ("FRINTL",FRINTL),
                ("FSINTL",FSINTL),
                ("FRINTS",FRINTS),
                ("FSINTS",FSINTS),
                ("CINTRL",CINTRL),
                ("CINTRS",CINTRS),
                ("CINTSL",CINTSL),
                ("CINTSS",CINTSS),
                ("MELFAC",MELFAC),
                ("CCFAC",CCFAC),
                ("LAIMLT",LAIMLT),
                ("SAIMLT",SAIMLT),
                ("GRDMLT",GRDMLT),
                ("MAXLQF",MAXLQF),
                ("KSNVP",KSNVP),
                ("SNODEN",SNODEN),
                ("GLMAX",GLMAX),
                ("CR",CR),
                ("GLMIN",GLMIN),
                ("RM",RM),
                ("R5",R5),
                ("CVPD",CVPD),
                ("TL",TL),
                ("T1",T1),
                ("T2",T2),
                ("TH",TH),
                ("MXKPL",MXKPL),
                ("MXRTLN",MXRTLN),
                ("inirlen",inirlen),
                ("inirdep",inirdep), # TODO(bernhard): only exported for derive_params_from_input_soil()
                ("rgrorate",rgrorate),
                ("rgroper",rgroper),
                ("FXYLEM",FXYLEM),
                ("PSICR",PSICR),
                ("RTRAD",RTRAD),
                ("NOOUTF",NOOUTF),
                ("FXYLEM",FXYLEM),
                ("inirlen",inirlen),
                ("NLAYER",NLAYER),
                ("nmat", nmat), # TODO(bernhard): only exported for derive_params_from_input_soil()
                ("ILAYER",ILAYER),
                ("QLAYER",QLAYER),
                ("IMODEL",IMODEL),
                ("RSSA",RSSA),
                ("RSSB",RSSB),
                ("INFEXP",INFEXP),
                ("BYPAR",BYPAR),
                ("QFPAR",QFPAR),
                ("QFFC",QFFC),
                ("IMPERV",IMPERV),
                ("DSLOPE",DSLOPE),
                ("LENGTH",LENGTH),
                ("DRAIN",DRAIN),
                ("GSC",GSC),
                ("GSP",GSP),
                ("DTIMAX",DTIMAX),
                ("DSWMAX",DSWMAX),
                ("DPSIMX",DPSIMX)])
end


end