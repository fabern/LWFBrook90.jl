# PLRF test

# parameters
NLAYER = 18
NOOUTF = true
p_PSICR = -1.071
p_RHOWG  = 0.00981

p_fu_DISPC = 15.309999999999999
p_fT_RXYLEM = 0.021794134560092007
p_fT_DAYLEN = 0.5882294500805856
p_STORAGEK = 2.4
p_CAPACITANCE = 2

p_fT_RROOTI = [0.12890668023781152, 0.8420843272927548, 1.7606151276361475, 0.2260820831279297, 2.365237126205569, 0.5164036888962881, 0.24262669526586808, 4.2686957705907815, 0.5003232600181334, 0.7363888076447798, 7.703990175011473, 0.5376625539004918, 13.903887230749579, 0.9342277043295377, 1.686060905291935, 3.702485894239593, 1.0e20, 1.0e20]
p_fT_ALPHA = [3.241223122887817e-6, 2.1972700449880832e-5, 4.649305474180868e-5, 6.202157469203119e-6, 6.741780128744168e-5, 1.504215510346679e-5, 7.543116539365436e-6, 0.00013970893129850516, 1.6946217416891348e-5, 2.6433992793479763e-5, 0.00028492081624199633, 2.098177478830277e-5, 0.0005736316442614217, 4.054597099172328e-5, 8.040709389401906e-5, 0.0002086824095753948, 1.0e20, 1.0e20]
p_fu_KK = [0.060319567894926644, 0.0018349724591501387, 0.0017919704261863665, 0.001594881398024343, 0.0010119540641645053, 0.0007580852302455766, 0.006439238222734642, 0.0061098201620713236, 0.0058971554955052794, 0.0057811748058775895, 0.005814662032754645, 0.005902350110427364, 0.0062841003374256785, 0.006754358757001102, 0.00866339546338736, 0.013252374258638084, 0.02153598439509493, 0.022984329029063416]

# states
u_PLPSI = -925.0426992896058
u_PLSTOR = 8.149914601420818
p_fu_PLPSIG = 0.5 * p_RHOWG * p_fu_DISPC * 1000 # plant gravity potential in kPa
u_aux_PLPSIT = u_PLPSI + p_fu_PLPSIG # total plant water potential

u_aux_PSITI = [-255.70471049220663, -253.60362659888597, -257.45419303925485, -276.94415358046166, -366.07789769269317, -436.9256859555962, -785.5409564592304, -813.2546924593253, -832.5180412248245, -844.1368100510414, -841.5056964126128, -834.4509381235944, -802.61590740218, -767.4672562807835, -656.81062974401, -504.22076841892334, -375.5185159956168, -363.2011845162336]     

p_fu_PTR = 2.820121590085702


## PLRFBYLAYER happens first in the day, then recompute PLPSI

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

    # top of loop for recalculation of root water uptake if more layers get flagged
    RI = zeros(NLAYER)
    RWU = zeros(NLAYER)
    while true
        # start loop with NEGFLAG = 0
        NEGFLAG = false

        ###
        SUM = 0
        for i = 1:NLAYER
            if (FLAG[i] == 0)
                RI[i] = p_fT_RROOTI[i] + p_fT_ALPHA[i] / p_fu_KK[i]
                SUM = SUM + 1.0 / RI[i]
            else
                RWU[i] = 0.0
            end
        end

        if SUM < 1E-20
            break
        end

        # distribute total plant refill to layers
        for  i = 1:NLAYER
            if FLAG[i] == 1
                RWU[i] = 0
            else
                RWU[i] = ((u_aux_PSITI[i] - u_aux_PLPSIT) / 1000) / (RI[i] + (1/p_STORAGEK))

                # check for any negative root water uptake
                if RWU[i] < -0.0000001
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

    u_aux_PLRF = sum((1 - p_fT_DAYLEN) .* RWU)

    aux_du_PLPSI = u_aux_PLRF / p_CAPACITANCE

    u_aux_PLPSIT = u_aux_PLPSIT + aux_du_PLPSI * 1000

    ##TBYLAYER

    FLAG = zeros(NLAYER+1) # add one for plant storage
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

    # Compute ATR, ATRANI and PLFL
    # top of loop for recalculation of transpiration if more layers get flagged
    RI = zeros(NLAYER)
    ATRANI=zeros(NLAYER)
    ATR = 0.0
    PLFL = 0.0
    while true
        # start loop with NEGFLAG = 0
        NEGFLAG = false

        ###
        SUM = 0
        for i = 1:NLAYER
            if (FLAG[i] == 0)
                RI[i] = p_fT_RROOTI[i] + p_fT_ALPHA[i] / p_fu_KK[i]
                SUM = SUM + 1.0 / RI[i]
                # 1 ./ sum((1 ./ (p_fT_RROOTI .+ p_fT_ALPHA ./ p_fu_KK))[FLAG .== 0])
            else
                ATRANI[i] = 0.0
            end
        end

        # check whether to include plant storage in transpiration
        if (FLAG[NLAYER+1] == 0)
            SUM = SUM + p_STORAGEK; # plant storage conductance
        else
            PLFL = 0.0
        end

        if SUM < 1E-20
                ATR = 0.
                PSIT = -1e10

                break
        else
                RT = 1.0 / SUM
        end

        ###
        # weighted mean soil water potential
        PSIT = 0
        for  i = 1:NLAYER
                if FLAG[i] == 0
                    PSIT = PSIT + RT * u_aux_PSITI[i] / RI[i]
                end
        end

        if (FLAG[NLAYER+1] == 0)
            PSIT = PSIT + RT * u_aux_PLPSIT * p_STORAGEK # add plant storage water potential weighted by storage conductance
        end

        # 1) compute available SUPPLY
        # soil water supply rate, assumed constant over day
        SUPPLY = (PSIT / 1000 - p_PSICR - p_RHOWG * p_fu_DISPC) / (RT + p_fT_RXYLEM)
                # PSIT                  in kPa
                # PSICR                 in MPa
                # RHOWG                 in MPa/m
                # DISPC                 in m
                # RT and RXYLEM must be in [MPa d/mm]  to get:
                # SUPPLY                in mm/day

        # 2) reduce actual transpiration rate ATR
        # transpiration rate limited by either PTR or SUPPLY
            # daytime, PTR is average of a half sine over daytime
            R = (2 / 3.14159) * (SUPPLY / p_fu_PTR)
            if R <= 0
                ATR = 0.
            elseif R < 1
                ATR = p_fu_PTR * (1 + R * acos(R) - sin(acos(R)))
            else
                ATR = p_fu_PTR
            end

        # distribute total transpiration rate to layers
        for  i = 1:NLAYER
            if FLAG[i] == 1
                ATRANI[i] = 0
            else
                ATRANI[i] = ((u_aux_PSITI[i]- PSIT) / 1000 + RT * ATR) / RI[i]

                # check for any negative transpiration losses
                if ATRANI[i] < -0.000001
                    NEGFLAG = true
                end
            end
        end
        
        # plant storage contribution to transpiration
        if FLAG[NLAYER+1] == 0
            PLFL = ((u_aux_PLPSIT - PSIT) / 1000 + RT * ATR) * p_STORAGEK # plant storage flux

            # prevent negative plant storage flux (downwards flux)
            if PLFL < -0.0000001
                NEGFLAG = true
            end
        else
            PLFL = 0.0
        end

        ###
        if NOOUTF && NEGFLAG
            # find layer with most negative transpiration and omit it
            IDEL = 0
            TRMIN = 0
            for i = 1:NLAYER
                if (ATRANI[i] < TRMIN)
                    TRMIN = ATRANI[i]
                    IDEL = i
                end
            end
            if (IDEL == 0) # then must be negative plant storage that is causing NEGFLAG
                IDEL = NLAYER + 1 # omit plant storage
            end
            FLAG[IDEL] = 1
        # repeat main loop with flagged layers excluded
        else
            break
        end
    end
