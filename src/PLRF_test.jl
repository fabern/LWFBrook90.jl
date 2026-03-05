# PLRF test

# parameters
NLAYER = 18
NOOUTF = true
p_PSICR = -1.2
p_RHOWG  = 0.00981

p_fu_DISPC = 15.309999999999999
p_fT_RXYLEM = 0.030110935023771795
p_STORAGEK = 2.4

p_fT_RROOTI = [0.18922725979118843, 1.5777545464328822, 3.15550909286576, 0.3506121214295292, 3.15550909286576, 0.6311018185731527, 0.29981487725250494, 5.140661116523625, 0.5140661116523632, 0.5711845685026258, 5.140661116523625, 0.5140306443290983, 1.0e20, 1.0e20, 1.0e20, 1.0e20, 1.0e20, 1.0e20]
p_fT_ALPHA = [4.449722975909061e-6, 4.081571690637576e-5, 8.16314338127514e-5, 9.070159312527945e-6, 8.16314338127514e-5, 1.6326286762550303e-5, 8.29114540322345e-6, 0.000148700916759153, 1.4870091675915321e-5, 1.652232408435036e-5, 0.000148700916759153, 1.6948256755129938e-5, 1.0e20, 1.0e20, 1.0e20, 1.0e20, 1.0e20, 1.0e20]
p_fu_KK = [0.8347500758890816, 1.7302441181614523, 1.7950816714979076, 2.0288716034793195, 2.214511326095429, 2.33435703026686, 2.9458413149201084, 3.110505411144997, 3.2303238396850866, 3.237788753868193, 3.271521794500662, 3.3481375476285025, 3.4618910896727257, 3.5865050270697925, 3.741260624437348, 3.864312570810623, 4.234818387090108, 6.374197755868146]

# states
u_PLPSI = -67.06848142139278
u_PLSTOR = 9.865863037157213

u_aux_PSITI = [-6.643350000000001, -7.0848, -7.23195, -7.72245, -8.212950000000001, -8.50725, -9.4392, -10.174949999999999, -10.714500000000005, -11.646449999999998, -12.136949999999999, -13.117949999999999, -14.098949999999999, -15.129000000000005, -17.091000000000005, -19.543499999999998, -22.977, -25.429499999999997]
u_aux_PSITI = [-10.965759640411276, -10.708957369674906, -10.63162137182339, -10.40969613540536, -10.421834751862578, -10.440034468700688, -10.659110847666074, -11.131535294617441, -11.493202601706562, -12.414422606119016, -12.856928365763736, -13.731830740692434, -14.562332417208012, -15.436412759527094, -17.216472102208936, -19.532749344196127, -22.595161923090423, -23.62216149806149]

    ## modified TBYLAYER for PLRF
    FLAG = zeros(NLAYER)
    for i = 1:NLAYER
        if p_fT_RROOTI[i] > 1E+15
            # This layer has no roots
            FLAG[i] = 1
        elseif (NOOUTF && u_aux_PSITI[i] / 1000. <= p_PSICR)
            # This layer has no outflow from roots to soil
            FLAG[i] = 1
        else
            # This layer has roots connected to soil
            FLAG[i]= 0
        end
    end

    # top of loop for recalculation of transpiration if more layers get flagged
    RI = zeros(NLAYER)
    RWU = zeros(NLAYER)

    PLRF = 0.0
    while true
        # start loop with NEGFLAG = 0
        NEGFLAG = false

        ###
        SUM = 0.0
        for i = 1:NLAYER
            if (FLAG[i] == 0)
                RI[i] = p_fT_RROOTI[i] + p_fT_ALPHA[i] / p_fu_KK[i]
                SUM = SUM + 1.0 / RI[i]
            else
                RWU[i] = 0.0
            end
        end

        if SUM < 1E-20
                PLRF = 0.0
                break
        end

        # distribute total plant refill to layers
        for  i = 1:NLAYER
            if FLAG[i] == 1
                RWU[i] = 0
            else
                RWU[i] = ((u_aux_PSITI[i] - u_PLPSI) / 1000 - 0.5 * p_RHOWG * p_fu_DISPC) / (RI[i] + (1/p_STORAGEK))

                # check for any negative root water uptake
                if RWU[i] < -0.000001
                    NEGFLAG = true
                end
            end
        end

        ###
        if NOOUTF && NEGFLAG
            # find layer with most negative root water uptake and omit it
            IDEL = 0
            RWUMIN = 0
            for i = 1:NLAYER
                if (RWU[i] < RWUMIN)
                    RWUMIN = RWU[i]
                    IDEL = i
                end
            end
            FLAG[IDEL] = 1
        # repeat main loop with flagged layers excluded
        else
            break
        end
    end
    
    sum(0.5 .* RWU)