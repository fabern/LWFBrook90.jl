using LWFBrook90
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

    @test LWFBrook90.PET.PM(150., -0.43, 0.06, 146., 19640.) ≈ 0.58512224
    @test LWFBrook90.PET.PM(150., 0.30, 0.06, 23.8, 578.) ≈ 14.0411894
end

@testset "KPT.KPT_SOILPAR_Mvg1d" begin
    p_soil1 = LWFBrook90.KPT.KPT_SOILPAR_Mvg1d(;
        p_THICK  = [40, 40, 120],
        p_STONEF = [0.010, 0.175, 0.175],
        p_THSAT  = [0.714, 0.668, 0.656],
        p_Kθfc   = [2.0, 2.0, 2.0],
        p_KSAT   = [24864, 12881, 10516],
        p_MvGα   = [1147,  1274,  1215],
        p_MvGn   = [1.051225, 1.051052, 1.051055],
        p_MvGl   = [4.6703, 4.4782, 4.5016],
        p_θr     = [0.069, 0.069, 0.069])

    p_soil2 = LWFBrook90.KPT.KPT_SOILPAR_Mvg1d(;
        p_THICK  = 3 .* [40, 40, 120],
        p_STONEF = [0.010, 0.175, 0.175],
        p_THSAT  = [0.714, 0.668, 0.656],
        p_Kθfc   = [2.0, 2.0, 2.0],
        p_KSAT   = [24864, 12881, 10516],
        p_MvGα   = [1147,  1274,  1215],
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
    @test p_soil1.p_SWATMX  ≈ [28.2744, 22.044, 64.944]
    @test p_soil1.p_WETF    ≈ [0.924422821232832, 0.9364352826531043, 0.9400828692513616]
    @test p_soil1.p_PsiCrit ≈ [-2.001039697148224e72, -7.555263344289398e71, -5.276582923029774e71]

    @test p_soil2.p_PSIG   ≈ [-0.5886, -1.7658, -4.1202000000000005]
    @test p_soil2.p_SWATMX ≈ [84.82319999999999, 66.132, 194.832]
    @test p_soil2.p_SWATMX ≈ p_soil1.p_SWATMX * 3
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

@testset "KPT.KPT_SOILPAR_Ch1d" begin
    # data from Ecoshift "b90v44data/SCsl.txt"
    p_soil1 = LWFBrook90.KPT.KPT_SOILPAR_Ch1d(;
        p_THICK  = [100,100],
        p_STONEF = [0,0],
        p_THSAT  = [0.435,0.435],
        p_PSIF   = [-7.90,-7.90],
        p_THETAF = [.266,.266],
        p_KF     = [5.50,5.50],
        p_BEXP   = [4.90,4.90],
        p_WETINF = [0.920, 0.29])

    p_soil2 = LWFBrook90.KPT.KPT_SOILPAR_Ch1d(;
        p_THICK  = 2 .* [100,100],
        p_STONEF = [0,0],
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
    @test p_soil1.p_SWATMX  ≈ [43.5, 43.5]
    @test p_soil1.p_WETF    ≈ [0.6114942528735633, 0.6114942528735633]
    @test p_soil1.p_PsiCrit ≈ [-4.7813122937946886e17, -4.7813122937946886e17]

    @test p_soil2.p_PSIG   ≈ [-0.9809999999999999, -2.9429999999999996]
    @test p_soil2.p_SWATMX ≈ [87.0, 87.0]
    @test p_soil2.p_SWATMX ≈ p_soil1.p_SWATMX * 2
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

@testset "run_example" begin
    example_result = LWFBrook90.run_example()

    idx_of_state_variables = 1:(7+(example_result["NLAYER"]-1))
    @test example_result["solution"].u[10][idx_of_state_variables] ≈
        [0.0, 0.0, 0.0, 0.829558, 0.0,
        0.041478, 23.641250, 16.321178,
        48.165395, 94.120618, 122.816595,
        150.917200, 1.914875]
        # [0.0, 0.0, 0.0, 0.8295578096607609, 0.0,
        # 0.041477890483038043, 23.641250144447312, 16.321178493182103,
        # 48.16539508679621, 94.12061823602788, 122.81659477698743,
        # 150.91719950931895, 1.9148754251864477]
    @test example_result["solution"].u[100][idx_of_state_variables] ≈
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 24.688027, 19.479705,
        58.713465, 123.462959, 145.403060,
        150.920137, 1.915018]
        # [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 24.688027405420044, 19.479705329409544,
        # 58.71346519060473, 123.46295863521023, 145.40306041260098,
        # 150.92013723079336, 1.9150179532438305]
    @test example_result["solution"].u[200][idx_of_state_variables] ≈
        [0.0, 1.387779e-17, 0.034885, 0.620356, 0.0,
        0.031018, 25.304292, 20.501754, 59.193256,
        115.708186, 155.423754, 176.530520, 1.9159532]
        # [0.0, 1.3877787807814457e-17, 0.034884656647269655, 0.6203560102512926, 0.0,
        # 0.03101780051256463, 25.304292411701322, 20.50175401792435, 59.19325594445119,
        # 115.70818579715453, 155.42375372071783, 176.53052033584265, 1.9159531711918132]
end