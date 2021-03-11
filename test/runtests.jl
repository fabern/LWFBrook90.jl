using LWFBrook90
using Test: @test, @testset

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

# @testset "Module PET" begin
@testset "run_example" begin
    example_result = LWFBrook90.run_example()
    idx_of_state_variables = 1:(7+(example_result["NLAYER"]-1))
    @test example_result["solution"].u[10][idx_of_state_variables] ≈
        [0.0, 0.0, 0.0, 0.8295578096607609, 0.0,
        0.041477890483038043, 23.641250144447312, 16.321178493182103,
        48.16539508679621, 94.12061823602788, 122.81659477698743,
        150.91719950931895, 1.9148754251864477]
    @test example_result["solution"].u[100][idx_of_state_variables] ≈
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 24.688027405420044, 19.479705329409544,
        58.71346519060473, 123.46295863521023, 145.40306041260098,
        150.92013723079336, 1.9150179532438305]
    @test example_result["solution"].u[200][idx_of_state_variables] ≈
        [0.0, 1.3877787807814457e-17, 0.034884656647269655, 0.6203560102512926, 0.0,
        0.03101780051256463, 25.304292411701322, 20.50175401792435, 59.19325594445119,
        115.70818579715453, 155.42375372071783, 176.53052033584265, 1.9159531711918132]
end
