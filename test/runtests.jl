using LWFBrook90
using DataFrames
using Test: @testset, @test, @test_throws

# @testset "Module WAT" begin
@testset "WAT.KKMEAN" begin
    @test LWFBrook90.WAT.KKMEAN(50.,  50., 1., 1.) ≈ 50
    @test LWFBrook90.WAT.KKMEAN(50.,  50., 1., 5.) ≈ 50
    @test LWFBrook90.WAT.KKMEAN(10., 100., 1., 1.) ≈ 31.6227766
    @test LWFBrook90.WAT.KKMEAN(10., 100., 1., 5.) ≈ 14.6779926
    @test LWFBrook90.WAT.KKMEAN(10., 100., 5., 1.) ≈ 68.1292069
end

# @testset "Module PET" begin
@testset "PET.ESAT" begin
    @test LWFBrook90.PET.ESAT(0)[1] ≈ 0.61078 #(0.61078, 0.044448927)
    @test LWFBrook90.PET.ESAT(0)[2] ≈ 0.044448927 #(0.61078, 0.044448927)
    @test LWFBrook90.PET.ESAT(25)[1] ≈ 3.1674898302368564  #(3.1674898302368564, 0.18866467946037985)
    @test LWFBrook90.PET.ESAT(25)[2] ≈ 0.18866467946037985 #(3.1674898302368564, 0.18866467946037985)

    # using Plots
    # plot(x -> x, x -> LWFBrook90.PET.ESAT(x)[1], -30, 45,
    #     xlabel = "T °C",
    #     ylabel = "ES [kPa] or dES/dT [kPa/°C]",
    #     label = "Saturation vapour pressure (ES)")
    # plot!(x -> x, x -> LWFBrook90.PET.ESAT(x)[2], -30, 45,
    #     label = "DELTA (dES/dT)")

    @test LWFBrook90.PET.PM(150., -0.43, 0.06, 146., 19640.) ≈ 0.58512224
    @test LWFBrook90.PET.PM(150., 0.30, 0.06, 23.8, 578.) ≈ 14.0411894
end

@testset "KPT.KPT_SOILPAR_Mvg1d" begin
    p_soil1 = KPT_SOILPAR_Mvg1d(;
        p_THICK  = [40., 40., 120.],
        p_STONEF = [0.010, 0.175, 0.175],
        p_THSAT  = [0.714, 0.668, 0.656],
        p_Kθfc   = [2.0, 2.0, 2.0],
        p_KSAT   = [24864., 12881., 10516.],
        p_MvGα   = [1147.,  1274.,  1215.],
        p_MvGn   = [1.051225, 1.051052, 1.051055],
        p_MvGl   = [4.6703, 4.4782, 4.5016],
        p_θr     = [0.069, 0.069, 0.069])

    p_soil2 = KPT_SOILPAR_Mvg1d(;
        p_THICK  = 3 .* [40., 40., 120.],
        p_STONEF = [0.010, 0.175, 0.175],
        p_THSAT  = [0.714, 0.668, 0.656],
        p_Kθfc   = [2.0, 2.0, 2.0],
        p_KSAT   = [24864., 12881., 10516.],
        p_MvGα   = [1147.,  1274.,  1215.],
        p_MvGn   = [1.051225, 1.051052, 1.051055],
        p_MvGl   = [4.6703, 4.4782, 4.5016],
        p_θr     = [0.069, 0.069, 0.069])

    @test p_soil1.p_THICK   ≈ [40, 40, 120]
    @test p_soil1.p_STONEF  ≈ [0.01, 0.175, 0.175]
    @test p_soil1.p_THSAT   ≈ [0.714, 0.668, 0.656]
    @test p_soil1.p_Kθfc    ≈ [2.0, 2.0, 2.0]
    @test p_soil1.p_KSAT    ≈ [24864, 12881, 10516]
    @test p_soil1.p_MvGα    ≈ [1147, 1274, 1215]
    @test p_soil1.p_MvGn    ≈ [1.051225, 1.051052, 1.051055]
    @test p_soil1.p_MvGl    ≈ [4.6703, 4.4782, 4.5016]
    @test p_soil1.p_θr      ≈ [0.069, 0.069, 0.069]
    @test p_soil1.p_PSIF    ≈ [-0.03210102700011719, -0.02096528992202748, -0.019805546964157438]
    @test p_soil1.p_THETAF  ≈ [0.6652527196951767, 0.6299247343092094, 0.6208286442505493]
    @test p_soil1.p_PSIG    ≈ [-0.19619999999999999, -0.5886, -1.3734]
    @test p_soil1.p_SWATMAX  ≈ [28.2744, 22.044, 64.944]
    @test p_soil1.p_WETF    ≈ [0.924422821232832, 0.9364352826531043, 0.9400828692513616]
    @test p_soil1.p_PsiCrit ≈ [-2.001039697148224e72, -7.555263344289398e71, -5.276582923029774e71]

    @test p_soil2.p_PSIG   ≈ [-0.5886, -1.7658, -4.1202000000000005]
    @test p_soil2.p_SWATMAX ≈ [84.82319999999999, 66.132, 194.832]
    @test p_soil2.p_SWATMAX ≈ p_soil1.p_SWATMAX * 3
    @test p_soil2.p_THICK  ≈ p_soil1.p_THICK * 3

    @test_throws AssertionError LWFBrook90.KPT.KPT_SOILPAR_Mvg1d(;
        p_THICK  = [1,1],
        p_STONEF = [1,1],
        p_THSAT  = [1,1],
        p_Kθfc   = [1,1],
        p_KSAT   = [1,1],
        p_MvGα   = [1,1],
        p_MvGn   = [1,1],
        p_MvGl   = [1,1],
        p_θr     = [1])
end

@testset "discretization" begin
    input_soil_discretization = DataFrame(Upper_m = -(0.:8.), Lower_m = -(1.:9.),
                               Rootden_ = (1:9) ./ 10, uAux_PSIM_init_kPa = -6.0,
                               u_delta18O_init_mUr = NaN, u_delta2H_init_mUr = NaN)
    input_soil_horizons = DataFrame(HorizonNr = 1:5, Upper_m = -[0.,3,8,10,15], Lower_m = -[3.,8,10,15,20],
                    ths_volFrac = (11:15) ./ 20, thr_volFrac = 0.069,
                    gravel_volFrac = 0.9:-0.1:0.5,
                    ksat_mmDay = 50, alpha_perMeter = 10:14, npar_ = 0, tort_ = 0
                    )
    IDEPTH = 45 # mm
    QDEPTH = 0  # mm
    INITRDEP = 10
    RGRORATE = 10
    FLAG_MualVanGen = 1
    # FLAG_MualVanGen = 0

    soil_disc = LWFBrook90.discretize_soil_params(
        input_soil_horizons, input_soil_discretization, [], IDEPTH, QDEPTH, INITRDEP, RGRORATE, FLAG_MualVanGen)

    @test soil_disc["NLAYER"] == 10
    @test soil_disc["THICK"]  ≈ [45.0, 955.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]
    @test soil_disc["STONEF"] ≈ [0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7]
    @test soil_disc["PSIM_init"] ≈ [-6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0]
    @test soil_disc["frelden"] ≈ [0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    @test soil_disc["QLAYER"] == 0
    @test soil_disc["ILAYER"] == 1
    #soil_disc["tini"] .= 0.0
    @test soil_disc["PAR"][!,"θs"] ≈ [0.55, 0.55, 0.55, 0.55, 0.6, 0.6, 0.6, 0.6, 0.6, 0.65]
end

@testset "KPT.KPT_SOILPAR_Ch1d" begin
    # data from Ecoshift "b90v44data/SCsl.txt"
    p_soil1 = KPT_SOILPAR_Ch1d(;
        p_THICK  = [100.,100.],
        p_STONEF = [0.,0.],
        p_THSAT  = [0.435,0.435],
        p_PSIF   = [-7.90,-7.90],
        p_THETAF = [.266,.266],
        p_KF     = [5.50,5.50],
        p_BEXP   = [4.90,4.90],
        p_WETINF = [0.920, 0.29])

    p_soil2 = KPT_SOILPAR_Ch1d(;
        p_THICK  = 2 .* [100.,100.],
        p_STONEF = [0.,0.],
        p_THSAT  = [0.435,0.435],
        p_PSIF   = [-7.90,-7.90],
        p_THETAF = [.266,.266],
        p_KF     = [5.50,5.50],
        p_BEXP   = [4.90,4.90],
        p_WETINF = [0.920, 0.29])

    @test p_soil1.p_THICK   ≈ [100, 100]
    @test p_soil1.p_STONEF  ≈ [0, 0]
    @test p_soil1.p_THSAT   ≈ [0.435, 0.435]
    @test p_soil1.p_PSIF    ≈ [-7.9, -7.9]
    @test p_soil1.p_THETAF  ≈ [0.266, 0.266]
    @test p_soil1.p_KF      ≈ [5.5, 5.5]
    @test p_soil1.p_BEXP    ≈ [4.9, 4.9]
    @test p_soil1.p_WETINF  ≈ [ 0.92, 0.29]
    @test p_soil1.p_CHM     ≈ [95.73164439624084, -6667.137100808987]
    @test p_soil1.p_CHN     ≈ [0.7806060606060607,0.3545656945751017]
    @test p_soil1.p_PSIG    ≈ [-0.49049999999999994, -1.4714999999999998]
    @test p_soil1.p_SWATMAX  ≈ [43.5, 43.5]
    @test p_soil1.p_WETF    ≈ [0.6114942528735633, 0.6114942528735633]
    @test p_soil1.p_PsiCrit ≈ [-4.7813122937946886e17, -4.7813122937946886e17]

    @test p_soil2.p_PSIG   ≈ [-0.9809999999999999, -2.9429999999999996]
    @test p_soil2.p_SWATMAX ≈ [87.0, 87.0]
    @test p_soil2.p_SWATMAX ≈ p_soil1.p_SWATMAX * 2
    @test p_soil2.p_THICK  ≈ p_soil1.p_THICK * 2


    @test_throws AssertionError LWFBrook90.KPT.KPT_SOILPAR_Ch1d(;
        p_THICK  = [1,1],
        p_STONEF = [1,1],
        p_THSAT  = [1,1],
        p_PSIF   = [1,1],
        p_THETAF = [1,1],
        p_KF     = [1,1],
        p_BEXP   = [1,1],
        p_WETINF = [1])
end

# TODO(bernhard): while below testset works consistently on MacBook Pro,
#                 in the GithubActions the values are slightly different.
#                 TODO: investigate why. (DiffEq.jl, integrative testing, set some seed?)

# @testset "run_example" begin
#     example_result = LWFBrook90.run_example()

#     idx_of_state_variables = 1:(7+(example_result["NLAYER"]-1))
#     @test example_result["solution"].u[10][idx_of_state_variables] ≈
#         [   0.0,
#             0.0,
#             0.0,
#             0.82955780,
#             0.0,
#             0.04147789,
#            23.68344511,
#            16.31494834,
#            48.15140732,
#            94.10625789,
#           122.80998540,
#           150.91516960,
#             1.91485613]
#         # [   0.0,
#         #     0.0,
#         #     0.0,
#         #     0.8295578096607609,
#         #     0.0,
#         #     0.041477890483038043,
#         #    23.683445116860266,
#         #    16.314948346294006,
#         #    48.15140732674836,
#         #    94.10625789421519,
#         #   122.80998540999502,
#         #   150.9151696005563,
#         #     1.9148561318569484]
#     @test example_result["solution"].u[100][idx_of_state_variables] ≈
#         [  0.0,
#            0.0,
#            0.0,
#            0.0,
#            0.0,
#            0.0,
#           25.14944841,
#           20.10942668,
#           58.82141224,
#          120.76261447,
#          149.55895926,
#          150.86467817,
#            1.91398950]
#         # [  0.0,
#         #    0.0,
#         #    0.0,
#         #    0.0,
#         #    0.0,
#         #    0.0,
#         #   25.14944841419801,
#         #   20.109426682979937,
#         #   58.8214122420572,
#         #  120.76261447765783,
#         #  149.55895926520608,
#         #  150.86467817418463,
#         #    1.9139895016342252]
#     @test example_result["solution"].u[200][idx_of_state_variables] ≈
#         [  0.0,
#            1.38777878e-17,
#            0.03488465,
#            0.62035601,
#            0.0,
#            0.03101780,
#           25.49412260,
#           20.33552067,
#           60.86255825,
#          113.55496111,
#          155.02076671,
#          181.33392684,
#            1.91378684]
#         # [  0.0,
#         #    1.3877787807814457e-17,
#         #    0.034884656647269655,
#         #    0.6203560102512926,
#         #    0.0,
#         #    0.03101780051256463,
#         #   25.494122609817502,
#         #   20.33552067169021,
#         #   60.86255825819409,
#         #  113.55496111973302,
#         #  155.02076671472267,
#         #  181.3339268423993,
#         #    1.913786847187571]
# end