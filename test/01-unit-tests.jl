# Unit tests

# - _Unit testing_ asserts that individual pieces of a project work as expected. (developers
#       perspective)
# - _Integration testing_ asserts that they fit together as expected. Also known as
#       _functional tests_, they cover entire use cases (user perspective). For LWFBrook90.jl
#       these are tests that are compared to e.g. LWFBrook90R or Hydrus.
# - _Regression testing_ asserts that behavior is unchanged over time. Also known as
#       _reference tests_.

# NOTE: locally, i.e. not on CI system, one might need to do manually cd("test")
if basename(pwd()) != "test"; cd("test"); end
@testset "Soil Hydraulics" begin
    # Test a single element of MualemVanGenuchtenSHP
    shp = LWFBrook90.MualemVanGenuchtenSHP(; p_THSAT = 1.0, p_θr = 1.0, p_MvGα = 1.0, p_MvGn = 1.0, p_KSAT = 1.0,
                                p_MvGl = 1.0, p_STONEF = 1.0)
    @test shp.p_THSAT ≈ 1.0
    @test shp.p_θr ≈ 1.0
    @test shp.p_MvGα ≈ 1.0
    @test shp.p_MvGn ≈ 1.0
    @test shp.p_KSAT ≈ 1.0
    @test shp.p_MvGl ≈ 1.0
    @test shp.p_STONEF ≈ 1.0

    # Test generating an array of MualemVanGenuchtenSHP'S
    soil_horizons = LWFBrook90.read_path_soil_horizons(
        "test-assets/DAV-2020/input-files/DAV_LW1_def_soil_horizons.csv");

    @test [shp.p_THSAT  for shp in soil_horizons.shp] ≈ [0.3786, 0.3786, 0.3786]
    @test [shp.p_θr     for shp in soil_horizons.shp] ≈ [0.0, 0.0, 0.0]
    @test [shp.p_MvGα   for shp in soil_horizons.shp] ≈ [20.387, 20.387, 20.387]
    @test [shp.p_MvGn   for shp in soil_horizons.shp] ≈ [1.2347, 1.2347, 1.2347]
    @test [shp.p_KSAT   for shp in soil_horizons.shp] ≈ [2854.91, 2854.91, 2854.91]
    @test [shp.p_MvGl   for shp in soil_horizons.shp] ≈ [-3.339, -3.339, -3.339]
    @test [shp.p_STONEF for shp in soil_horizons.shp] ≈ [0.175, 0.375, 0.75]
end

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
    # @test p_soil1.p_Kθfc    ≈ [2.0, 2.0, 2.0]
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

    @test_logs (:warn, r"\[a,b\] is not a bracketing interval") min_level = Logging.Warn match_mode=:any KPT_SOILPAR_Mvg1d(;
        p_THICK  = [40., 40., 120.],
        p_STONEF = [0.010, 0.175, 0.175],
        p_THSAT  = [0.714, 0.668, 0.656],
        # p_Kθfc   = [2.0, 2.0, 2.0],
        p_Kθfc   = [1000002.0, 1000002.0, 1000002.0],
        p_KSAT   = [24864., 12881., 10516.],
        p_MvGα   = [1147.,  1274.,  1215.],
        p_MvGn   = [1.051225, 1.051052, 1.051055],
        p_MvGl   = [4.6703, 4.4782, 4.5016],
        p_θr     = [0.069, 0.069, 0.069])
end

@testset "discretization" begin
    input_soil_discretization = DataFrame(Upper_m = -(0.:8.), Lower_m = -(1.:9.),
                                Rootden_ = (1:9) ./ 10, uAux_PSIM_init_kPa = -6.0,
                                u_delta18O_init_permil = NaN, u_delta2H_init_permil = NaN)
    # input_soil_horizons = LWFBrook90.read_path_soil_horizons(
    #     "test-assets/DAV-2020/input-files/DAV_LW1_def_soil_horizons.csv");
    input_soil_horizons =
        DataFrame(HorizonNr = 1:5, Upper_m = -[0.,3,8,10,15], Lower_m = -[3.,8,10,15,20],
                    shp = LWFBrook90.MualemVanGenuchtenSHP(;
                        p_THSAT = 0.55, p_θr = 0.069, p_MvGα = 12, p_MvGn = 0, p_KSAT = 50,
                        p_MvGl = 0, p_STONEF = 0.7))
                        # p_THSAT = (11:15) ./ 20, p_θr = 0.069, p_MvGα = 10:14, p_MvGn = 0, p_KSAT = 50,
                        # p_MvGl = 0, p_STONEF = 0.9:-0.1:0.5))
    IDEPTH_m = 0.045 # m
    QDEPTH_m = 0.0  # m
    INITRDEP = 10
    RGRORATE = 10
    FLAG_MualVanGen = 1
    # FLAG_MualVanGen = 0

    soil_disc = LWFBrook90.refine_soil_discretization(
        input_soil_horizons, input_soil_discretization, [], IDEPTH_m, QDEPTH_m)

    @test soil_disc["NLAYER"] == 10
    @test soil_disc["THICK"]  ≈ [45.0, 955.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]
    @test soil_disc["PSIM_init"] ≈ [-6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0]
    keys(soil_disc)
    @test soil_disc["final_Rootden_"] ≈ [0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    @test soil_disc["QLAYER"] == 0
    @test soil_disc["ILAYER"] == 1
    @test [shp.p_STONEF for shp in soil_disc["SHP"]] ≈ [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7]
    # @test [shp.p_STONEF for shp in soil_disc["SHP"]] ≈ [0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7]
    @test [shp.p_THSAT for shp in soil_disc["SHP"]] ≈ [0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55]
    # @test [shp.p_THSAT for shp in soil_disc["SHP"]] ≈ [0.55, 0.55, 0.55, 0.55, 0.6, 0.6, 0.6, 0.6, 0.6, 0.65]
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

# @testset "adding-soil-layers" begin
Δz_m_data = [
            [fill(0.5, 4);],
            [fill(0.02, 100);],
            [fill(0.01, 200);]
        ]
@testset "adding-soil-layers (Δz_m = $(first(Δz_m)))" for Δz_m in Δz_m_data # source: https://stackoverflow.com/a/63871951
    continuous_SPAC = SPAC("test-assets/Hammel-2001/input-files-ISO/",
                      "Hammel_loam-NLayer-103-RESET=FALSE";
                      simulate_isotopes = false);
    f1 = (Δz_m) -> LWFBrook90.Rootden_beta_(0.97, Δz_m = Δz_m)  # function for root density as f(Δz)
    f2 = (Δz_m) -> fill(-6.3, length(Δz_m))          # function for initial conditions as f(Δz)

    soil_discretization = discretize_soil(;
        Δz_m = Δz_m,
        Rootden_ = f1,
        uAux_PSIM_init_kPa = f2,
        u_delta18O_init_permil = ifelse.(cumsum(Δz_m) .<= 0.2, -13., -10.),
        u_delta2H_init_permil  = ifelse.(cumsum(Δz_m) .<= 0.2, -95., -70.))

    ####################
    ## Discretize soil parameters and interpolate discretized root distribution
    # Define refinement of grid with soil_output_depths
    if length(Δz_m) == 4
        soil_output_depths = [-0.91, -1.01, -1.11, -1.21, -1.31]
    else
        soil_output_depths = soil_output_depths = zeros(Float64, 0)
    end
    soil_discr =
        LWFBrook90.refine_soil_discretization(
            continuous_SPAC.soil_horizons,
            soil_discretization,
            soil_output_depths,
            continuous_SPAC.params[:IDEPTH_m],
            continuous_SPAC.params[:QDEPTH_m])

    # Interpolate discretized root distribution in time
    p_fT_RELDEN = LWFBrook90.HammelKennel_transient_root_density(;
        timepoints = continuous_SPAC.meteo_forcing.p_days,
        p_AGE      = continuous_SPAC.canopy_evolution.p_AGE,
        p_INITRDEP = continuous_SPAC.params[:INITRDEP],
        p_INITRLEN = continuous_SPAC.params[:INITRLEN],
        p_RGROPER_y  = continuous_SPAC.params[:RGROPER],
        p_RGRORATE_m_per_y = continuous_SPAC.params[:RGRORATE],
        p_THICK         = soil_discr["THICK"],
        final_Rootden_profile = soil_discr["final_Rootden_"]);
    ####################

    ####################
    # Define parameters for differential equation
    p = define_LWFB90_p(continuous_SPAC, soil_discr, p_fT_RELDEN);
    # using Plots
    # hline([0; cumsum(p.p_THICK)], yflip = true, xticks = false,
    #     title = "N_layer = "*string(p.NLAYER))
   ####################

    ####################
    # Define state vector u for DiffEq.jl
    # a) allocation of u0
    u0 = define_LWFB90_u0(;simulate_isotopes = continuous_SPAC.solver_options.simulate_isotopes,
                          compute_intermediate_quantities = continuous_SPAC.solver_options.compute_intermediate_quantities,
                          NLAYER = soil_discr["NLAYER"]);
    ####################

    # Check if defined layers correspond to requested
    @test p.p_soil.NLAYER == length(Δz_m) + 2*length(soil_output_depths) + 1 # +1 because we needed to add one at 0.005 m for the IDEPTH_m
end

# TODO(bernhard): include unit tests of specific functions, e.g. during development of the Hammel-2001 infiltration test
#                 it was noticed, that ISVP (snow evaporation rate could become negative, therby generating intercepted snow out of nowwhere...)
# using LWFBrook90.EVP
# # INTER24(, , , , , p_FRINTS, , , , , MONTHN)
# @run LWFBrook90.EVP.INTER24(0.0, # p_fT_RFAL or p_fT_SFAL
#     -0.8811, # p_fu_PINT
#     1.7520, # p_fu_LAI
#     1.0, # p_fu_SAI
#     0.0, # p_FRINTL or p_FSINTL
#     0.0, # p_FRINTS or p_FSINTS
#     0.6, # p_CINTRL or p_CINTSL
#     0.6, # p_CINTRS or p_CINTSS
#     4.0, # p_DURATN
#     0.0, # u_INTR or u_INTS
#     1  # MONTHN
# ) # returns (aux_du_RINT, aux_du_IRVP)

# using LWFBrook90.EVP
# LWFBrook90.EVP.INTER24
# # INTER24(, , , , , p_FRINTS, , , , , MONTHN)
# @run LWFBrook90.EVP.INTER24(0.0, # p_fT_RFAL or p_fT_SFAL
#     -0.8811, # p_fu_PINT
#     1.7520, # p_fu_LAI
#     1.0, # p_fu_SAI
#     0.0, # p_FRINTL or p_FSINTL
#     0.0, # p_FRINTS or p_FSINTS
#     0.6, # p_CINTRL or p_CINTSL
#     0.6, # p_CINTRS or p_CINTSS
#     4.0, # p_DURATN
#     0.0, # u_INTR or u_INTS
#     1  # MONTHN
# ) # returns (aux_du_RINT, aux_du_IRVP)