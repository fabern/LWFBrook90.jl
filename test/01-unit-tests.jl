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
    p_soil1 = LWFBrook90.KPT.KPT_SOILPAR_Mvg1d(;
        p_THICK  = [40., 40., 120.],
        p_STONEF = [0.010, 0.175, 0.175],
        p_THSAT  = [0.714, 0.668, 0.656],
        p_Kθfc   = [2.0, 2.0, 2.0],
        p_KSAT   = [24864., 12881., 10516.],
        p_MvGα   = [1147.,  1274.,  1215.],
        p_MvGn   = [1.051225, 1.051052, 1.051055],
        p_MvGl   = [4.6703, 4.4782, 4.5016],
        p_θr     = [0.069, 0.069, 0.069])

    p_soil2 = LWFBrook90.KPT.KPT_SOILPAR_Mvg1d(;
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

    @test_logs (:warn, r"\[a,b\] is not a bracketing interval") min_level = Logging.Warn match_mode=:any LWFBrook90.KPT.KPT_SOILPAR_Mvg1d(;
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

@testset "reading_forcing" begin
    path_meteoveg = "test-assets/BEA-2016/input-files/BEA2016-reset-FALSE_meteoveg.csv"
    meteo_forcing, reference_date = LWFBrook90.read_path_meteoveg(path_meteoveg)
    @test all(diff(meteo_forcing.days) .== 1.0)

    path_meteoveg = "test-assets/BEA-2016/input-files/BEA2016-reset-FALSE_meteoveg_withGaps.csv"
    @test_throws AssertionError meteo_forcing, reference_date = LWFBrook90.read_path_meteoveg(path_meteoveg)
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
    p_soil1 = LWFBrook90.KPT.KPT_SOILPAR_Ch1d(;
        p_THICK  = [100.,100.],
        p_STONEF = [0.,0.],
        p_THSAT  = [0.435,0.435],
        p_PSIF   = [-7.90,-7.90],
        p_THETAF = [.266,.266],
        p_KF     = [5.50,5.50],
        p_BEXP   = [4.90,4.90],
        p_WETINF = [0.920, 0.29])

    p_soil2 = LWFBrook90.KPT.KPT_SOILPAR_Ch1d(;
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
    f2 = (Δz_m) -> fill(-6.3, length(Δz_m))                     # function for initial conditions as f(Δz)
    f3 = (Δz_m) -> ifelse.(cumsum(Δz_m) .<= 0.2, -13., -10.)    # function for initial conditions as f(Δz)
    f4 = (Δz_m) -> ifelse.(cumsum(Δz_m) .<= 0.2, -95., -70.)    # function for initial conditions as f(Δz)

    soil_discretization = LWFBrook90.discretize_soil(;
        Δz_m = Δz_m,
        Rootden_ = f1,
        uAux_PSIM_init_kPa = f2,
        u_delta18O_init_permil = f3,
        u_delta2H_init_permil  = f4)

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
    p = LWFBrook90.define_LWFB90_p(continuous_SPAC, soil_discr, p_fT_RELDEN);
    # using Plots
    # hline([0; cumsum(p.p_THICK)], yflip = true, xticks = false,
    #     title = "N_layer = "*string(p.NLAYER))
   ####################

    ####################
    # Define state vector u for DiffEq.jl
    # a) allocation of u0
    u0 = LWFBrook90.define_LWFB90_u0(;simulate_isotopes = continuous_SPAC.solver_options.simulate_isotopes,
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

@testset "root-model (Δz_m = $(first(Δz_m)))" for Δz_m in Δz_m_data # source: https://stackoverflow.com/a/63871951
    continuous_SPAC = SPAC("test-assets/Hammel-2001/input-files-ISO",
                      "Hammel_loam-NLayer-27-RESET=FALSE";
                      simulate_isotopes = false);

    continuous_SPAC.root_distribution = (beta = 0.95, )

    simulation = discretize(continuous_SPAC, Δz = (thickness_m = Δz_m,
                                      functions = (rootden = nothing,  # function for root density as f(Δz)
                                                    PSIM_init = (Δz_m) -> fill(-6.3, length(Δz_m)),                      # function for initial conditions as f(Δz)
                                                    δ18Ο_init = (Δz_m) -> ifelse.(cumsum(Δz_m) .<= 0.2, -13., -10.),    # function for initial conditions as f(Δz)
                                                    δ2Η_init = (Δz_m) -> ifelse.(cumsum(Δz_m) .<= 0.2, -95., -70.))));  # function for initial conditions as f(Δz))))

    # plot(simulation.soil_discretization.Rootden_, simulation.soil_discretization.Lower_m)
    if ([0.5, 0.5, 0.5, 0.5] == Δz_m)
        @test simulation.soil_discretization.Rootden_ ≈ [0.018461100494465737, 0.018461100494465737, 0.0014204889211275828, 0.00010929948491700703, 8.410046164697427e-6]
    elseif ([fill(0.02, 100);] == Δz_m)
        @test simulation.soil_discretization.Rootden_ ≈ [0.048750000000000016, 0.048750000000000016, 0.043996875000000046, 0.039707179687500094, 0.03583572966796872, 0.03234174602534179, 0.02918842578787101, 0.026342554273553587, 0.023774155231882088, 0.021456175096773555, 0.01936419802483813, 0.017476188717416408, 0.015772260317468367, 0.014234464936515145, 0.012846604605204925, 0.011594060656197502, 0.0104636397422182, 0.009443434867351885, 0.008522699967785152, 0.007691736720926101, 0.006941792390635748, 0.006264967632548801, 0.005654133288375274, 0.005102855292758668, 0.004605326901714724, 0.0041563075287975315, 0.00375106754473975, 0.0033853384591276403, 0.003055267959362673, 0.002757379333324872, 0.002488534848325674, 0.0022459027006138665, 0.002026927187304073, 0.0018293017865418926, 0.0016509448623540646, 0.0014899777382745838, 0.0013447049087927376, 0.0012135961801855166, 0.0010952705526173778, 0.0009884816737372182, 0.0008921047105478475, 0.0008051245012694053, 0.0007266248623956084, 0.0006557789383120904, 0.0005918404918266451, 0.0005341360438735343, 0.00048205777959586804, 0.00043505714608527146, 0.00039263907434194945, 0.0003543567645936663, 0.0003198069800456893, 0.0002886257994912933, 0.00026048478404089304, 0.00023508751759687696, 0.00021216648463123766, 0.00019148025237963884, 0.00017281092777265972, 0.0001559618623148129, 0.00014075558073911587, 0.00012703191161700378, 0.00011464630023438893, 0.0001034682859615832, 9.338012808024487e-5, 8.427556559248428e-5, 7.605869794719e-5, 6.864297489733717e-5, 6.195028484484721e-5, 5.591013207245643e-5, 5.045889419541538e-5, 4.553915201138681e-5, 4.1099084690243703e-5, 3.709192393291927e-5, 3.3475461349519176e-5, 3.0211603867902337e-5, 2.7265972490808643e-5, 2.460754017291622e-5, 2.2208305006099494e-5, 2.0042995267977037e-5, 1.8088803229343586e-5, 1.6325144914508538e-5, 1.4733443285341874e-5, 1.329693256502118e-5, 1.2000481639906635e-5, 1.083043468003142e-5, 9.77446729871767e-6, 8.821456737084787e-6, 7.96136470526676e-6, 7.185131646492149e-6, 6.484581310917115e-6, 5.8523346331273984e-6, 5.281732006368056e-6, 4.766763135821694e-6, 4.30200373002787e-6, 3.882558366363753e-6, 3.5040089256255236e-6, 3.1623680553649614e-6, 2.85403716998589e-6, 2.5757685458982493e-6, 2.324631112715636e-6, 2.097979579174236e-6, 1.893426570254153e-6]
    elseif ([fill(0.01, 200);] == Δz_m)
        @test simulation.soil_discretization.Rootden_ ≈ [0.050000000000000044, 0.050000000000000044, 0.04749999999999999, 0.04512500000000008, 0.04286875000000001, 0.0407253125, 0.03868904687500008, 0.03675459453125007, 0.034916864804687475, 0.033171021564453174, 0.031512470486230404, 0.029936846961919006, 0.02844000461382301, 0.027018004383131844, 0.02566710416397522, 0.024383748955776552, 0.023164561507987735, 0.022006333432588177, 0.020906016760958934, 0.019860715922910943, 0.01886768012676543, 0.01792429612042712, 0.017028081314405807, 0.016176677248685434, 0.01536784338625119, 0.01459945121693862, 0.013869478656091783, 0.01317600472328695, 0.012517204487122902, 0.011891344262766612, 0.011296777049628282, 0.010731938197146795, 0.010195341287289605, 0.009685574222925042, 0.00920129551177884, 0.00874123073618982, 0.008304169199380373, 0.007888960739411255, 0.0074945127024408364, 0.00711978706731875, 0.006763797713952857, 0.0064256078282550755, 0.006104327436842416, 0.005799111065000306, 0.005509155511750241, 0.005233697736162779, 0.004972012849354557, 0.00472341220688699, 0.004487241596542457, 0.004262879516715445, 0.004049735540879618, 0.0038472487638356867, 0.0036548863256438135, 0.003472142009361745, 0.0032985349088936466, 0.0031336081634488755, 0.0029769277552764706, 0.0028280813675126693, 0.0026866772991369636, 0.0025523434341802043, 0.002424726262471144, 0.0023034899493475924, 0.0021883154518802517, 0.0020788996792862058, 0.001974954695321829, 0.001876206960555793, 0.0017823966125279922, 0.0016932767819016759, 0.0016086129428065643, 0.0015281822956662028, 0.0014517731808828538, 0.0013791845218387166, 0.0013102252957468696, 0.0012447140309593818, 0.0011824783294115404, 0.001123354412940869, 0.0010671866922938866, 0.0010138273576791867, 0.0009631359897952496, 0.000914979190305476, 0.0008692302307902189, 0.0008257687192506635, 0.000784480283288147, 0.0007452562691236952, 0.0007079934556675216, 0.000672593782884201, 0.0006389640937399799, 0.0006070158890529864, 0.0005766650946003038, 0.0005478318398702831, 0.0005204402478767856, 0.0004944182354829074, 0.00046969732370882866, 0.00044621245752340943, 0.0004239018346471335, 0.0004027067429148712, 0.0003825714057690277, 0.00036344283548073175, 0.0003452706937066008, 0.00032800715902125965, 0.00031160680107023, 0.00029602646101667407, 0.0002812251379658015, 0.00026716388106762246, 0.0002538056870141636, 0.0002411154026634721, 0.00022905963253028183, 0.00021760665090386766, 0.00020672631835860766, 0.00019639000244064952, 0.00018657050231862815, 0.0001772419772027023, 0.00016837987834261714, 0.0001599608844253808, 0.000151962840204245, 0.00014436469819389952, 0.00013714646328433222, 0.00013028914012003234, 0.00012377468311397521, 0.00011758594895838748, 0.00011170665151039039, 0.00010612131893494858, 0.00010081525298821781, 9.5774490338707e-5, 9.098576582178275e-5, 8.643647753070471e-5, 8.211465365426385e-5, 7.800892097142853e-5, 7.410847492295147e-5, 7.04030511767817e-5, 6.688289861789265e-5, 6.353875368703132e-5, 6.0361816002663105e-5, 5.734372520249664e-5, 5.447653894241622e-5, 5.175271199531206e-5, 4.91650763955187e-5, 4.670682257579273e-5, 4.437148144698089e-5, 4.2152907374637394e-5, 4.004526200585001e-5, 3.804299890552976e-5, 3.614084896030878e-5, 3.433380651229889e-5, 3.261711618673946e-5, 3.0986260377341424e-5, 2.943694735846325e-5, 2.7965099990590048e-5, 2.656684499102724e-5, 2.523850274149808e-5, 2.397657760433436e-5, 2.27777487241676e-5, 2.1638861288031386e-5, 2.0556918223557652e-5, 1.9529072312396423e-5, 1.85526186967655e-5, 1.7624987761921673e-5, 1.6743738373903305e-5, 1.590655145511377e-5, 1.5111223882424696e-5, 1.4355662688259052e-5, 1.3637879553884957e-5, 1.2955985576157403e-5, 1.2308186297271817e-5, 1.1692776982541453e-5, 1.1108138133320011e-5, 1.0552731226742829e-5, 1.0025094665411238e-5, 9.523839932024103e-6, 9.047647935522818e-6, 8.595265538646757e-6, 8.165502261792135e-6, 7.757227148741386e-6, 7.369365791265459e-6, 7.000897501718839e-6, 6.650852626521875e-6, 6.3183099953123545e-6, 6.002394495552288e-6, 5.702274770702509e-6, 5.417161032195139e-6, 5.146302980540973e-6, 4.888987831574987e-6, 4.644538440068402e-6, 4.412311517931755e-6, 4.191695942123985e-6, 3.98211114494007e-6, 3.7830055877874358e-6, 3.593855308348104e-6, 3.414162542902943e-6, 3.243454415713387e-6, 3.0812816950165356e-6, 2.9272176103045666e-6, 2.7808567296672138e-6, 2.6418138933115287e-6, 2.50972319848497e-6, 2.384237038688397e-6, 2.265025186742875e-6, 2.1517739273724246e-6, 2.0441852309760478e-6, 1.9419759694772054e-6, 1.8448771710311007e-6]
    end
end

@testset "root-model (Δz_m = $(first(Δz_m)))" for Δz_m in Δz_m_data # source: https://stackoverflow.com/a/63871951
    continuous_SPAC = SPAC("test-assets/Hammel-2001/input-files-ISO",
                      "Hammel_loam-NLayer-27-RESET=FALSE";
                      simulate_isotopes = false);

    continuous_SPAC.root_distribution = (beta = 0.95, maxRootDepth_m = 0.5)

    simulation = discretize(continuous_SPAC, Δz = (thickness_m = Δz_m,
                                      functions = (rootden = nothing,  # function for root density as f(Δz)
                                                    PSIM_init = (Δz_m) -> fill(-6.3, length(Δz_m)),                      # function for initial conditions as f(Δz)
                                                    δ18Ο_init = (Δz_m) -> ifelse.(cumsum(Δz_m) .<= 0.2, -13., -10.),    # function for initial conditions as f(Δz)
                                                    δ2Η_init = (Δz_m) -> ifelse.(cumsum(Δz_m) .<= 0.2, -95., -70.))));  # function for initial conditions as f(Δz))))

    # plot(simulation.soil_discretization.Rootden_, simulation.soil_discretization.Lower_m)
    if ([0.5, 0.5, 0.5, 0.5] == Δz_m)
        @test simulation.soil_discretization.Rootden_ ≈ [0.018461100494465737, 0.018461100494465737, 0.0, 0.0, 0.0]
    elseif ([fill(0.02, 100);] == Δz_m)
        @test simulation.soil_discretization.Rootden_ ≈ [0.048750000000000016, 0.048750000000000016, 0.043996875000000046, 0.039707179687500094, 0.03583572966796872, 0.03234174602534179, 0.02918842578787101, 0.026342554273553587, 0.023774155231882088, 0.021456175096773555, 0.01936419802483813, 0.017476188717416408, 0.015772260317468367, 0.014234464936515145, 0.012846604605204925, 0.011594060656197502, 0.0104636397422182, 0.009443434867351885, 0.008522699967785152, 0.007691736720926101, 0.006941792390635748, 0.006264967632548801, 0.005654133288375274, 0.005102855292758668, 0.004605326901714724, 0.004156307528797476, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    elseif ([fill(0.01, 200);] == Δz_m)
        @test simulation.soil_discretization.Rootden_ ≈ [0.050000000000000044, 0.050000000000000044, 0.04749999999999999, 0.04512500000000008, 0.04286875000000001, 0.0407253125, 0.03868904687500008, 0.03675459453125007, 0.034916864804687475, 0.033171021564453174, 0.031512470486230404, 0.029936846961919006, 0.02844000461382301, 0.027018004383131844, 0.02566710416397522, 0.024383748955776552, 0.023164561507987735, 0.022006333432588177, 0.020906016760958934, 0.019860715922910943, 0.01886768012676543, 0.01792429612042712, 0.017028081314405807, 0.016176677248685434, 0.01536784338625119, 0.01459945121693862, 0.013869478656091783, 0.01317600472328695, 0.012517204487122902, 0.011891344262766612, 0.011296777049628282, 0.010731938197146795, 0.010195341287289605, 0.009685574222925042, 0.00920129551177884, 0.00874123073618982, 0.008304169199380373, 0.007888960739411255, 0.0074945127024408364, 0.00711978706731875, 0.006763797713952857, 0.0064256078282550755, 0.006104327436842416, 0.005799111065000306, 0.005509155511750241, 0.005233697736162779, 0.004972012849354557, 0.00472341220688699, 0.004487241596542457, 0.004262879516715445, 0.004049735540879507, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    end
end
