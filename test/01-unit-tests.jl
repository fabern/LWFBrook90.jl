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
    @test shp.p_MvGm ≈ 1.0 - 1.0 / 1.0
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
    @test [shp.p_MvGm   for shp in soil_horizons.shp] ≈ [0.19008666072730207, 0.19008666072730207, 0.19008666072730207]
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
        p_MvGm   = [0.048728863944445866, 0.04857228757473475, 0.04857500321105945],
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
        p_MvGm   = [0.048728863944445866, 0.04857228757473475, 0.04857500321105945],
        p_MvGl   = [4.6703, 4.4782, 4.5016],
        p_θr     = [0.069, 0.069, 0.069])

    @test p_soil1.p_THICK   ≈ [40, 40, 120]
    @test p_soil1.p_STONEF  ≈ [0.01, 0.175, 0.175]
    @test p_soil1.p_THSAT   ≈ [0.714, 0.668, 0.656]
    # @test p_soil1.p_Kθfc    ≈ [2.0, 2.0, 2.0]
    @test p_soil1.p_KSAT    ≈ [24864, 12881, 10516]
    @test p_soil1.p_MvGα    ≈ [1147, 1274, 1215]
    @test p_soil1.p_MvGn    ≈ [1.051225, 1.051052, 1.051055]
    @test p_soil1.p_MvGm    ≈ [0.048728863944445866, 0.04857228757473475, 0.04857500321105945]
    @test p_soil1.p_MvGl    ≈ [4.6703, 4.4782, 4.5016]
    @test p_soil1.p_θr      ≈ [0.069, 0.069, 0.069]
    @test p_soil1.p_PSIF    ≈ [-0.03210102700011719, -0.02096528992202748, -0.019805546964157438]
    @test p_soil1.p_THETAF  ≈ [0.6652527196951767, 0.6299247343092094, 0.6208286442505493]
    @test p_soil1.p_PSIG    ≈ [-0.19619999999999999, -0.5886, -1.3734]
    @test p_soil1.p_SWATMAX  ≈ [28.2744, 22.044, 64.944]
    @test p_soil1.p_WETF    ≈ [0.924422821232832, 0.9364352826531043, 0.9400828692513616]

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
        p_MvGm   = [0,0],
        p_MvGl   = [1,1],
        p_θr     = [1])

    @test_throws r"Computed invalid p_WETF in KPT_SOILPAR_Mvg1d" LWFBrook90.KPT.KPT_SOILPAR_Mvg1d(;
        p_THICK  = [40., 40., 120.],
        p_STONEF = [0.010, 0.175, 0.175],
        p_THSAT  = [0.714, 0.668, 0.656],
        # p_Kθfc   = [2.0, 2.0, 2.0],
        p_Kθfc   = [1000002.0, 1000002.0, 1000002.0],
        p_KSAT   = [24864., 12881., 10516.],
        p_MvGα   = [1147.,  1274.,  1215.],
        p_MvGn   = [1.051225, 1.051052, 1.051055],
        p_MvGm   = [0.048728863944445866, 0.04857228757473475, 0.04857500321105945],
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
    refined_soil_disc, IDEPTH_idx, QDEPTH_idx =
        LWFBrook90.refine_soil_discretization(
            # modifiedSPAC.soil_discretization.Δz,
            input_soil_discretization,
            input_soil_horizons,
            [],
            IDEPTH_m, QDEPTH_m)
    final_soil_discr = LWFBrook90.map_soil_horizons_to_discretization(input_soil_horizons, refined_soil_disc)

    @test nrow(final_soil_discr) == 10
    THICK_m        = final_soil_discr[!,"Upper_m"] - final_soil_discr[!,"Lower_m"] # thickness of soil layer [m]
    THICK          = 1000*(THICK_m)                                 # thickness of soil layer [mm]
    @test THICK  ≈ [45.0, 955.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]
    @test final_soil_discr[!,"uAux_PSIM_init_kPa"] ≈ [-6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0, -6.0]
    @test final_soil_discr[!,"Rootden_"] ≈ [0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    @test QDEPTH_idx == 0
    @test IDEPTH_idx == 1
    @test [shp.p_STONEF for shp in final_soil_discr.shp] ≈ [0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7]
    # @test [shp.p_STONEF for shp in final_soil_discr.shp] ≈ [0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7]
    @test [shp.p_THSAT for shp in final_soil_discr.shp] ≈ [0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55]
    # @test [shp.p_THSAT for shp in final_soil_discr.shp] ≈ [0.55, 0.55, 0.55, 0.55, 0.6, 0.6, 0.6, 0.6, 0.6, 0.65]
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

Δz_m_data = [
            [fill(0.5, 4);],
            [fill(0.02, 100);],
            [fill(0.01, 200);]
        ]
@testset "adding-soil-layers (Δz_m = $(first(Δz_m)))" for Δz_m in Δz_m_data # source: https://stackoverflow.com/a/63871951
    parametrizedSPAC = loadSPAC("test-assets/Hammel-2001/input-files-ISO/",
                      "Hammel_loam-NLayer-103-RESET=FALSE";
                      simulate_isotopes = false);
    # TODO: Note just because of this test we need the fct discretize_soil()
    #       and we set ε to 0.005m i.e. 5mm. Probably 10mm would be close enough.
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
    # Define refinement of grid with soil_output_depths_m
    if length(Δz_m) == 4
        soil_output_depths_m = [-0.95, -1.05, -1.15, -1.25, -1.35]
    else
        soil_output_depths_m = soil_output_depths_m = zeros(Float64, 0)
    end
    refined_soil_discretizationDF, IDEPTH_idx, QDEPTH_idx =
        LWFBrook90.refine_soil_discretization(
            # parametrizedSPAC.soil_discretization.df,
            soil_discretization,
            parametrizedSPAC.pars.soil_horizons,
            soil_output_depths_m,
            parametrizedSPAC.pars.params[:IDEPTH_m],
            parametrizedSPAC.pars.params[:QDEPTH_m],
            ε = 0.005)

    # Discretize the model in space as `soil_discretization`
    final_soil_discretizationDF = LWFBrook90.map_soil_horizons_to_discretization(parametrizedSPAC.pars.soil_horizons, refined_soil_discretizationDF)#computational_grid)

    # Update soil_discretization in underlying SPAC model
    parametrizedSPAC.soil_discretization = (
        Δz = final_soil_discretizationDF.Upper_m - final_soil_discretizationDF.Lower_m,
        df = final_soil_discretizationDF)

    ## c) Derive time evolution of aboveground vegetation based on parameter from SPAC
    canopy_evolution_relative = LWFBrook90.generate_canopy_timeseries_relative(
        parametrizedSPAC.pars.canopy_evolution,
        days = parametrizedSPAC.forcing.meteo["p_days"],
        reference_date = parametrizedSPAC.reference_date)
    canopy_evolutionDF = LWFBrook90.make_absolute_from_relative(
                aboveground_relative          = canopy_evolution_relative,
                p_MAXLAI                      = parametrizedSPAC.pars.params[:MAXLAI],
                p_SAI_baseline_               = parametrizedSPAC.pars.params[:SAI_baseline_],
                p_DENSEF_baseline_            = parametrizedSPAC.pars.params[:DENSEF_baseline_],
                p_AGE_baseline_yrs            = parametrizedSPAC.pars.params[:AGE_baseline_yrs],
                p_HEIGHT_baseline_m           = parametrizedSPAC.pars.params[:HEIGHT_baseline_m])
    ####################

    ####################
    ## d) Interpolate vegetation parameter in time for use as parameters
    # Aboveground: LAI, SAI, DENSEF, HEIGHT, AGE
    vegetation_fT = LWFBrook90.interpolate_aboveground_veg(canopy_evolutionDF.AboveGround)

    ## Interpolate discretized root distribution in time
        # b) Make root growth module on final discretized soil...
    vegetation_fT["p_RELDEN"] = LWFBrook90.HammelKennel_transient_root_density(;
        timepoints         = parametrizedSPAC.forcing.meteo["p_days"],
        AGE_at_timepoints  = vegetation_fT["p_AGE"].(parametrizedSPAC.forcing.meteo["p_days"]),
        p_INITRDEP         = parametrizedSPAC.pars.params[:INITRDEP],
        p_INITRLEN         = parametrizedSPAC.pars.params[:INITRLEN],
        p_RGROPER_y        = parametrizedSPAC.pars.params[:RGROPER],
        p_RGRORATE_m_per_y = parametrizedSPAC.pars.params[:RGRORATE],
        p_THICK               = 1000*parametrizedSPAC.soil_discretization.Δz,
        final_Rootden_profile = parametrizedSPAC.soil_discretization.df.Rootden_)
    ####################

    ####################
    # Define parameters for differential equation
    p = LWFBrook90.define_LWFB90_p(parametrizedSPAC, vegetation_fT, IDEPTH_idx, QDEPTH_idx);

    # using Plots
    # hline([0; cumsum(p.p_soil.p_THICK)], yflip = true, xticks = false,
    #     title = "N_layer = "*string(p.NLAYER))
    ####################

    ####################
    # Define state vector u for DiffEq.jl and initial states u0
        # state vector: GWAT,INTS,INTR,SNOW,CC,SNOWLQ,SWATI
    # a) allocation of u0
    u0 = LWFBrook90.define_LWFB90_u0(;simulate_isotopes = parametrizedSPAC.solver_options.simulate_isotopes,
                          compute_intermediate_quantities = parametrizedSPAC.solver_options.compute_intermediate_quantities,
                          NLAYER = nrow(parametrizedSPAC.soil_discretization.df));
    # b) initialization of u0
    LWFBrook90.init_LWFB90_u0!(;u0=u0, parametrizedSPAC=parametrizedSPAC, p_soil=p.p_soil);
    ####################

    # Check if defined layers correspond to requested
    expected_NLAYER = length(Δz_m) + 2*length(soil_output_depths_m) + 1 # +1 because we needed to add one at 0.005 m for the IDEPTH_m
    @test nrow(refined_soil_discretizationDF) == expected_NLAYER
    @test p.p_soil.NLAYER                     == expected_NLAYER
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

@testset "root-model beta (Δz_m = $(first(Δz_m)))" for Δz_m in Δz_m_data # source: https://stackoverflow.com/a/63871951
    parametrizedSPAC = loadSPAC("test-assets/Hammel-2001/input-files-ISO",
                      "Hammel_loam-NLayer-27-RESET=FALSE";
                      simulate_isotopes = false,
                      Δz_thickness_m = Δz_m,
                      root_distribution = (beta = 0.95, ),
                      IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -9.0, delta2H_init_permil = -11.0))
    # parametrizedSPAC.pars.root_distribution
    simulation = setup(parametrizedSPAC, ε = 0.005)

    # plot(simulation.parametrizedSPAC.soil_discretization.df.Rootden_, simulation.parametrizedSPAC.soil_discretization.df.Lower_m)
    if ([0.5, 0.5, 0.5, 0.5] == Δz_m)
        # @test simulation.parametrizedSPAC.soil_discretization.df.Rootden_ ≈ [0.018461100494465737, 0.018461100494465737, 0.0014204889211275828, 0.00010929948491700703, 8.410046164697427e-6]
        @test simulation.parametrizedSPAC.soil_discretization.df.Rootden_ ≈ [0.7201999283091463, 0.25792436725268547, 0.020201681877684945, 0.0015544179126264956, 0.00011960464785673769]
    elseif ([fill(0.02, 100);] == Δz_m)
        @test simulation.parametrizedSPAC.soil_discretization.df.Rootden_[1:10:100] ≈ [0.09207591746299604, 0.03520806629944542, 0.01262159612357775, 0.004524664528628884, 0.0016220285371347674, 0.0005814743963076271, 0.00020845038531732468, 7.472652865688285e-5, 2.6788408553974765e-5, 9.603267350270948e-6]
    elseif ([fill(0.01, 200);] == Δz_m)
        @test simulation.parametrizedSPAC.soil_discretization.df.Rootden_[1:10:200] ≈ [0.04823125873564055, 0.030012878588235002, 0.017969819063652963, 0.010759194464839009, 0.006441927162548269, 0.0038570197521006896, 0.002309340200954726, 0.0013826872835797437, 0.0008278659520944154, 0.0004956739262566578, 0.0002967782894672004, 0.00017769212466797057, 0.00010639083885043477, 6.370012521630651e-5, 3.8139618001095856e-5, 2.2835598145652136e-5, 1.367251613943244e-5, 8.186240464953907e-6, 4.901404559846817e-6, 2.93465196418323e-6]
    end
end

@testset "root-model beta cropped (Δz_m = $(first(Δz_m)))" for Δz_m in Δz_m_data # source: https://stackoverflow.com/a/63871951
    parametrizedSPAC = loadSPAC("test-assets/Hammel-2001/input-files-ISO",
                      "Hammel_loam-NLayer-27-RESET=FALSE";
                      simulate_isotopes = false,
                      Δz_thickness_m = Δz_m,
                      root_distribution = (beta = 0.95, z_rootMax_m = -0.5),
                      IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -9.0, delta2H_init_permil = -11.0))
    simulation = setup(parametrizedSPAC, ε = 0.005)

    # plot(simulation.soil_discretization.Rootden_, simulation.soil_discretization.Lower_m)
    if ([0.5, 0.5, 0.5, 0.5] == Δz_m)
        @test simulation.parametrizedSPAC.soil_discretization.df.Rootden_ ≈ [0.73630716625382, 0.26369283374617997, 0.0, 0.0, 0.0]
    elseif ([fill(0.02, 100);] == Δz_m)
        @test simulation.parametrizedSPAC.soil_discretization.df.Rootden_[1:5:50] ≈ [0.09899773754146445, 0.06322449003427368, 0.0378548376480284, 0.02266508962874636, 0.013570426391879014, 0.0081251155620333, 0.0, 0.0, 0.0, 0.0]
    elseif ([fill(0.01, 200);] == Δz_m)
        @test simulation.parametrizedSPAC.soil_discretization.df.Rootden_[1:10:100] ≈ [0.05204344024903263, 0.032385085823067206, 0.01939014716267554, 0.011609597363562048, 0.006951094791249151, 0.0041618772196682025, 0.0, 0.0, 0.0, 0.0]
    end
end
@testset "bare-minimum provided to loadSPAC" begin
    Δz_m = fill(0.1, 11)
    parametrizedSPAC = loadSPAC(
        # "examples/DAV2020-bare-minimum/", "DAV2020-minimal";
        "../examples/DAV2020-bare-minimum/", "DAV2020-minimal";
        simulate_isotopes = true,
        Δz_thickness_m = Δz_m,
        root_distribution = (beta = 0.77, z_rootMax_m = -0.5),
        IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -10.11111, delta2H_init_permil = -91.1111),
        canopy_evolution = (DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel    = 100,
                                        LAI_rel = (DOY_Bstart = 120,
                                            Bduration  = 21,
                                            DOY_Cstart = 270,
                                            Cduration  = 60,
                                            LAI_perc_BtoC = 95,
                                            LAI_perc_CtoB = 70)),
        storm_durations_h = [5.44, 5.44, 5.44, 5.44, 5.44, 5.44, 5.44, 5.44, 5.44, 5.44, 5.44, 5.44],
        IC_scalar = (amount = (u_GWAT_init_mm = 1.,
                               u_INTS_init_mm = 13.7,
                               u_INTR_init_mm = 0.,
                               u_SNOW_init_mm = 22.222,
                               u_CC_init_MJ_per_m2 = 0.101010,
                               u_SNOWLQ_init_mm =  0.),
                    d18O    = (u_GWAT_init_permil = -11.111,
                               u_INTS_init_permil = -12.222,
                               u_INTR_init_permil = -13.333,
                               u_SNOW_init_permil = -14.444),
                    d2H     = (u_GWAT_init_permil = -95.111,
                               u_INTS_init_permil = -95.222,
                               u_INTR_init_permil = -95.333,
                               u_SNOW_init_permil = -95.444)));

    @test_throws r"tspan \([0-9., ]*\) goes beyond input forcing data" setup(parametrizedSPAC, requested_tspan = (0., 400)) # forcing only defined from 0 to 364
    @test_throws AssertionError simulation = setup(parametrizedSPAC, soil_output_depths_m = [-1.0755, -1.096])
    simulation                             = setup(parametrizedSPAC, soil_output_depths_m = [-1.0755, -1.096], ε = 0.005);
    # Test soil discretization
    ## Δz
    actual_interfaces = cumsum(simulation.ODEProblem.p.p_soil.p_THICK)
    # we specified IDEPTH at 33.3cm, which means a layer should have been added to have an interface at 333cm:
    # further we specified a soil horizon ending at 3cm this requires another layer at 30mm
    # lastly we requested an output value at a depths of -1.0755, -1.094 m,
    #         hence for -1075.5mm two interfaces added (-1075.5mm, -1080.5mm) to have 5mm thick layer at that depth
    #         hence for -1094.0mm one interface added (-1094.0mm) and not (-1099.0mm) as we already had one at 1100mm which is close enought
    @test parametrizedSPAC.pars.params.IDEPTH_m ≈ 0.333

    @test sum(simulation.ODEProblem.p.p_soil.p_THICK[simulation.ODEProblem.p.p_INFRAC .!= 0.0]) ≈ 333.0

    @test any(actual_interfaces .- 0.00001 .<  333.0 .< actual_interfaces .+ 0.00001)
    @test any(actual_interfaces .- 0.00001 .<   30.0 .< actual_interfaces .+ 0.00001)
    @test any(actual_interfaces .- 0.00001 .< 1075.5 .< actual_interfaces .+ 0.00001)
    @test any(actual_interfaces .- 0.00001 .< 1080.5 .< actual_interfaces .+ 0.00001)
    @test any(actual_interfaces .- 0.00001 .< 1096.0 .< actual_interfaces .+ 0.00001)

    ## root_distribution
    simulation.parametrizedSPAC.soil_discretization.df.Rootden_ ≈
        [
         0.0926733195274138,
         0.0926733195274138,
         0.0067898780051124374,
         0.0004974726659130013,
         3.64482326699056e-5,
         3.64482326699056e-5,
         2.6704455456272315e-6,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0,
         0.0
        ]

    ## initial conditions
    (_, u_aux_PSIM, _, _, _) = # (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        LWFBrook90.KPT.derive_auxiliary_SOILVAR(simulation.ODEProblem.u0.SWATI.mm,
                                                simulation.ODEProblem.p.p_soil);
    @test all(u_aux_PSIM .≈ -7.0)
    @test all(simulation.ODEProblem.u0.SWATI.d18O .== fill(-10.11111, simulation.ODEProblem.p.p_soil.NLAYER))
    @test all(simulation.ODEProblem.u0.SWATI.d2H  .== fill(-91.1111, simulation.ODEProblem.p.p_soil.NLAYER))
    # Test canopy evolution
    @test simulation.ODEProblem.p.p_LAI(1) ≈ parametrizedSPAC.pars.params.MAXLAI * 70/100
    @test simulation.ODEProblem.p.p_LAI(180) ≈ parametrizedSPAC.pars.params.MAXLAI * 95/100
    # Test storm durations
    @test simulation.ODEProblem.p.p_DURATN == fill(5.44, 12)

    # Test scalar initial conditions
    @test parametrizedSPAC.pars.IC_scalar.u_INTS_init_mm[1] == 13.7
    @test simulation.ODEProblem.u0.INTS.mm                  == 13.7
    @test parametrizedSPAC.pars.IC_scalar.u_SNOW_init_mm[1] == 22.222
    @test simulation.ODEProblem.u0.SNOW.mm                  == 22.222
    @test parametrizedSPAC.pars.IC_scalar.u_CC_init_MJ_per_m2[1] == 0.101010
    @test simulation.ODEProblem.u0.CC.MJm2                       == 0.101010

    @test parametrizedSPAC.pars.IC_scalar.u_INTS_init_mm[2] == -12.222
    @test simulation.ODEProblem.u0.INTS.d18O                == -12.222
    @test parametrizedSPAC.pars.IC_scalar.u_SNOW_init_mm[2] == -14.444
    @test simulation.ODEProblem.u0.SNOW.d18O                == -14.444

    @test simulation.ODEProblem.u0.GWAT.d18O == -11.111
    @test simulation.ODEProblem.u0.INTS.d18O == -12.222
    @test simulation.ODEProblem.u0.INTR.d18O == -13.333
    @test simulation.ODEProblem.u0.SNOW.d18O == -14.444
    @test simulation.ODEProblem.u0.GWAT.d2H  == -95.111
    @test simulation.ODEProblem.u0.INTS.d2H  == -95.222
    @test simulation.ODEProblem.u0.INTR.d2H  == -95.333
    @test simulation.ODEProblem.u0.SNOW.d2H  == -95.444
end

@testset "remake-SPAC" begin
    Δz_m = fill(0.1, 11)
    parametrizedSPAC = loadSPAC(
            "../examples/DAV2020-full/", "DAV2020-full"; simulate_isotopes = true,
            # "../examples/DAV2020-bare-minimum/", "DAV2020-minimal"; simulate_isotopes = true,
            Δz_thickness_m = Δz_m,
            root_distribution = (beta = 0.77, z_rootMax_m = -0.5),
            IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -10.11111, delta2H_init_permil = -91.1111),
            canopy_evolution = (DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel    = 100,
                                            LAI_rel = (DOY_Bstart = 120,
                                                Bduration  = 21,
                                                DOY_Cstart = 270,
                                                Cduration  = 60,
                                                LAI_perc_BtoC = 95,
                                                LAI_perc_CtoB = 70)))#,
            # storm_durations_h = [5.44, 5.44, 5.44, 5.44, 5.44, 5.44, 5.44, 5.44, 5.44, 5.44, 5.44, 5.44],
            # IC_scalar = (amount = (u_GWAT_init_mm = 1.,
            #                        u_INTS_init_mm = 13.7,
            #                        u_INTR_init_mm = 0.,
            #                        u_SNOW_init_mm = 22.222,
            #                        u_CC_init_MJ_per_m2 = 0.101010,
            #                        u_SNOWLQ_init_mm =  0.),
            #             d18O    = (u_GWAT_init_permil = -11.111,
            #                        u_INTS_init_permil = -12.222,
            #                        u_INTR_init_permil = -13.333,
            #                        u_SNOW_init_permil = -14.444),
            #             d2H     = (u_GWAT_init_permil = -95.111,
            #                        u_INTS_init_permil = -95.222,
            #                        u_INTR_init_permil = -95.333,
            #                        u_SNOW_init_permil = -95.444)));

    # Discretized the baseline
    discrSPAC = setup(parametrizedSPAC);

    # Remake (modifies and (re)-discretizes)

    # TEST CHANGES TO INITIAL CONDITIONS ################################################################
    ICscalar_tochange = copy(discrSPAC.parametrizedSPAC.pars.IC_scalar)
    ICscalar_tochange.u_INTR_init_mm = [1.3, -10, -90.]
    ICscalar_tochange.u_GWAT_init_mm = [10, -10, -90.]
    remSPAC_0  = remakeSPAC(discrSPAC,
        IC_soil = (PSIM_init_kPa = -3.0, delta18O_init_permil = -15.55, ),
        IC_scalar = ICscalar_tochange);
    # test parametrizedSPAC:
    @test remSPAC_0.parametrizedSPAC.pars.IC_scalar.u_INTR_init_mm == [1.3, -10, -90.]
    @test remSPAC_0.parametrizedSPAC.pars.IC_scalar.u_GWAT_init_mm == [10, -10, -90.]
    @test remSPAC_0.parametrizedSPAC.pars.IC_soil.PSIM_init_kPa == -3.0
    @test remSPAC_0.parametrizedSPAC.pars.IC_soil.delta18O_init_permil == -15.55
    # test ODEProblem:
    @test remSPAC_0.ODEProblem.u0.INTR.mm    == 1.3
    @test remSPAC_0.ODEProblem.u0.INTR.d18O  == -10
    @test remSPAC_0.ODEProblem.u0.INTR.d2H   == -90
    @test remSPAC_0.ODEProblem.u0.GWAT.mm    == 10
    @test remSPAC_0.ODEProblem.u0.GWAT.d18O  == -10
    @test remSPAC_0.ODEProblem.u0.GWAT.d2H   == -90
    (u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
        LWFBrook90.KPT.derive_auxiliary_SOILVAR(
            remSPAC_0.ODEProblem.u0.SWATI.mm,
            remSPAC_0.ODEProblem.p.p_soil);
    @test all(u_aux_PSIM .≈ -3.0)
    @test all(remSPAC_0.ODEProblem.u0.SWATI.d18O .≈ -15.55)

    # TEST CHANGES TO SOIL HYDRAULICS ################################################################
    # all horizons proportionally (with `soil_horizons=(ths_ = 0.4)`)
        remSPAC_1  = remakeSPAC(discrSPAC, soil_horizons = (ths_ = 0.4,))
        # test parametrizedSPAC:
        @test names(discrSPAC.parametrizedSPAC.soil_discretization.df) == names(remSPAC_1.parametrizedSPAC.soil_discretization.df)
        @test remSPAC_1.parametrizedSPAC.pars.soil_horizons.shp[1].p_THSAT != discrSPAC.parametrizedSPAC.pars.soil_horizons.shp[1].p_THSAT
        @test remSPAC_1.parametrizedSPAC.pars.soil_horizons.shp[1].p_THSAT == 0.4
        # test ODEProblem:
        @test remSPAC_1.ODEProblem.p.p_soil.p_THSAT[1] == 0.4
        @test remSPAC_1.ODEProblem.p.p_soil.p_THSAT[end] != discrSPAC.parametrizedSPAC.pars.soil_horizons.shp[end].p_THSAT

        remSPAC_2  = remakeSPAC(discrSPAC, soil_horizons = (Ksat_mmday = 3854.9,))
        # test parametrizedSPAC:
        @test remSPAC_2.parametrizedSPAC.pars.soil_horizons.shp[1].p_KSAT != discrSPAC.parametrizedSPAC.pars.soil_horizons.shp[1].p_KSAT
        @test remSPAC_2.parametrizedSPAC.pars.soil_horizons.shp[1].p_KSAT == 3854.9
        # test ODEProblem:
        @test remSPAC_2.ODEProblem.p.p_soil.p_KSAT[1] == 3854.9
        @test remSPAC_2.ODEProblem.p.p_soil.p_KSAT[end] != discrSPAC.parametrizedSPAC.pars.soil_horizons.shp[end].p_KSAT

        remSPAC_3  = remakeSPAC(discrSPAC, soil_horizons = (alpha_per_m = 7.11,)) # we modify alpha as this scales h (it seems we are off by some orders in SCH)
        # test parametrizedSPAC:
        @test remSPAC_3.parametrizedSPAC.pars.soil_horizons.shp[1].p_MvGα != discrSPAC.parametrizedSPAC.pars.soil_horizons.shp[1].p_MvGα
        @test remSPAC_3.parametrizedSPAC.pars.soil_horizons.shp[1].p_MvGα == 7.11
        # test ODEProblem:
        @test remSPAC_3.ODEProblem.p.p_soil.p_MvGα[1] == 7.11
        @test remSPAC_3.ODEProblem.p.p_soil.p_MvGα[end] != discrSPAC.parametrizedSPAC.pars.soil_horizons.shp[end].p_MvGα
    # all horizons independently (with `soil_horizons=(ths_ = [0.4, 0.3, 0.3, 0.2])`) containing vectors for each soil_horizons
        @test_throws AssertionError remSPAC_1b  = remakeSPAC(discrSPAC, soil_horizons = (ths_ = [0.4, 0.3, 0.3, 0.2], ))
        remSPAC_1b  = remakeSPAC(discrSPAC, soil_horizons = (ths_ = [0.4, 0.3, 0.2], Ksat_mmday = 3801, ))
        # test parametrizedSPAC:
        @test names(discrSPAC.parametrizedSPAC.soil_discretization.df) == names(remSPAC_1b.parametrizedSPAC.soil_discretization.df)
        @test remSPAC_1b.parametrizedSPAC.pars.soil_horizons.shp[1].p_THSAT == 0.4
        @test remSPAC_1b.parametrizedSPAC.pars.soil_horizons.shp[2].p_THSAT == 0.3
        @test remSPAC_1b.parametrizedSPAC.pars.soil_horizons.shp[3].p_THSAT == 0.2

        @test remSPAC_1b.parametrizedSPAC.pars.soil_horizons.shp[1].p_KSAT ≈ 3801
        @test remSPAC_1b.parametrizedSPAC.pars.soil_horizons.shp[2].p_KSAT ≈ 3801
        @test remSPAC_1b.parametrizedSPAC.pars.soil_horizons.shp[3].p_KSAT ≈ 3801

        # test ODEProblem:
        @test all(remSPAC_1b.ODEProblem.p.p_soil.p_THSAT[[1]]     .== 0.4)
        @test all(remSPAC_1b.ODEProblem.p.p_soil.p_THSAT[[2:5; ]] .== 0.3)
        @test all(remSPAC_1b.ODEProblem.p.p_soil.p_THSAT[[6:12;]] .== 0.2)
        @test all(remSPAC_1b.ODEProblem.p.p_soil.p_KSAT .≈ 3801)

        @test_throws AssertionError remSPAC_3b  = remakeSPAC(discrSPAC, soil_horizons = (npar_ = [0.4, 0.3, 0.3, 0.2], ))
        remSPAC_3b  = remakeSPAC(discrSPAC, soil_horizons = (npar_ = [1.1, 1.2, 1.3], Ksat_mmday = 3801, ))
        # test parametrizedSPAC:
        @test names(discrSPAC.parametrizedSPAC.soil_discretization.df) == names(remSPAC_3b.parametrizedSPAC.soil_discretization.df)
        @test remSPAC_3b.parametrizedSPAC.pars.soil_horizons.shp[1].p_MvGn == 1.1
        @test remSPAC_3b.parametrizedSPAC.pars.soil_horizons.shp[2].p_MvGn == 1.2
        @test remSPAC_3b.parametrizedSPAC.pars.soil_horizons.shp[3].p_MvGn == 1.3

        @test remSPAC_3b.parametrizedSPAC.pars.soil_horizons.shp[1].p_MvGm == 1 - 1/1.1
        @test remSPAC_3b.parametrizedSPAC.pars.soil_horizons.shp[2].p_MvGm == 1 - 1/1.2
        @test remSPAC_3b.parametrizedSPAC.pars.soil_horizons.shp[3].p_MvGm == 1 - 1/1.3

    # TEST CHANGES TO FLOW ################################################################
    to_change = (DRAIN=.33, BYPAR=1, IDEPTH_m=0.67, INFEXP=0.33)
    remSPAC_4  = remakeSPAC(discrSPAC, params = to_change)
    # test parametrizedSPAC:
    @test remSPAC_4.parametrizedSPAC.pars.params[[:DRAIN, :BYPAR, :IDEPTH_m, :INFEXP]] != discrSPAC.parametrizedSPAC.pars.params[[:DRAIN, :BYPAR, :IDEPTH_m, :INFEXP]]
    @test remSPAC_4.parametrizedSPAC.pars.params[[:DRAIN, :BYPAR, :IDEPTH_m, :INFEXP]] == to_change
    # test ODEProblem:
    @test remSPAC_4.ODEProblem.p.p_DRAIN == to_change.DRAIN
    @test remSPAC_4.ODEProblem.p.p_BYPAR == to_change.BYPAR
    @test remSPAC_4.ODEProblem.p.p_INFRAC == LWFBrook90.WAT.INFPAR(
        to_change.INFEXP,
        sum(-to_change.IDEPTH_m .<= remSPAC_4.parametrizedSPAC.soil_discretization.df.Lower_m), # new ILAYER
        remSPAC_4.ODEProblem.p.p_soil)

    # TEST CHANGES TO POT TRANSPIRATION ################################################################
    to_change = (GLMAX=.00801, R5=235., CVPD=1.9,) # discrSPAC.parametrizedSPAC.pars.params[[:GLMAX, :R5, :CVPD]]
    remSPAC_5  = remakeSPAC(discrSPAC, params = to_change)
    # test parametrizedSPAC:
    @test all(values(remSPAC_5.parametrizedSPAC.pars.params[keys(to_change)]) .≈ values(Dict(pairs(to_change),)))
    # test ODEProblem:
    @test remSPAC_5.ODEProblem.p.p_GLMAX ≈ to_change.GLMAX
    @test remSPAC_5.ODEProblem.p.p_R5    ≈ to_change.R5
    @test remSPAC_5.ODEProblem.p.p_CVPD  ≈ to_change.CVPD

    # TEST CHANGES TO INTERCEPTION ################################################################
    to_change = (CINTRL=0.18,) # discrSPAC.parametrizedSPAC.pars.params[[:CINTRL]]
    remSPAC_6  = remakeSPAC(discrSPAC, params = to_change)
    # test parametrizedSPAC:
    @test all(values(remSPAC_6.parametrizedSPAC.pars.params[keys(to_change)]) .≈ values(Dict(pairs(to_change),)))
    # test ODEProblem:
    @test all(values(remSPAC_6.ODEProblem.p[Symbol.("p_".*String.(keys(to_change)))]) .≈ values(to_change))

    # TEST CHANGES TO ENERGY BALANCE ################################################################
    to_change = (ALB=0.15, ALBSN=0.7,) # discrSPAC.parametrizedSPAC.pars.params[[:ALB, :ALBSN]]
    remSPAC_7  = remakeSPAC(discrSPAC, params = to_change)
    # test parametrizedSPAC:
    @test remSPAC_7.parametrizedSPAC.pars.params[[:ALB, :ALBSN]] == to_change
    # test ODEProblem:
    @test all(values(remSPAC_7.ODEProblem.p[Symbol.("p_".*String.(keys(to_change)))]) .≈ values(to_change))

    # TEST CHANGES TO SOIL EVAPORATION
    to_change = (RSSA=720.,) # discrSPAC.parametrizedSPAC.pars.params[[:RSSA]]
    remSPAC_8  = remakeSPAC(discrSPAC, params = to_change)
    # test parametrizedSPAC:
    @test all(values(remSPAC_8.parametrizedSPAC.pars.params[keys(to_change)]) .≈ values(Dict(pairs(to_change),)))
    # test ODEProblem:
    @test all(values(remSPAC_8.ODEProblem.p[Symbol.("p_".*String.(keys(to_change)))]) .≈ values(to_change))

    # TEST CHANGES TO PLANT WATER SUPPLY ################################################################
    to_change = (PSICR=-1.6, FXYLEM=0.4, MXKPL=16.5) # discrSPAC.parametrizedSPAC.pars.params[[:PSICR, :FXYLEM, :MXKPL]]
    remSPAC_9  = remakeSPAC(discrSPAC, params = to_change)
    # test parametrizedSPAC:
    @test remSPAC_9.parametrizedSPAC.pars.params[[:PSICR, :FXYLEM, :MXKPL]] == to_change
    # test ODEProblem:
    @test all(values(remSPAC_9.ODEProblem.p[Symbol.("p_".*String.(keys(to_change)))]) .≈ values(to_change))

    # TEST CHANGES TO PLANT ################################################################
    to_change = (MAXLAI=9.999,) # discrSPAC.parametrizedSPAC.pars.params[[:MAXLAI]]
    remSPAC_10  = remakeSPAC(discrSPAC, params = to_change)
    # test parametrizedSPAC:
    @test all(values(remSPAC_10.parametrizedSPAC.pars.params[keys(to_change)]) .≈ values(Dict(pairs(to_change),)))
    # test ODEProblem:
    if discrSPAC.parametrizedSPAC.pars.canopy_evolution isa DataFrame
        max_relative_LAI = maximum(discrSPAC.parametrizedSPAC.pars.canopy_evolution.LAI_rel)/100
    elseif discrSPAC.parametrizedSPAC.pars.canopy_evolution isa NamedTuple
        max_relative_LAI = discrSPAC.parametrizedSPAC.pars.canopy_evolution.LAI_rel.LAI_perc_BtoC/100
    else
        error("...")
    end
    # @test maximum(remSPAC_10.ODEProblem.p.p_LAI.(1:365)) ≈ discrSPAC.parametrizedSPAC.pars.params[:MAXLAI] * max_relative_LAI
    @test maximum(remSPAC_10.ODEProblem.p.p_LAI.(1:365)) ≈ to_change.MAXLAI * max_relative_LAI
    # plot(remSPAC_10.ODEProblem.p.p_LAI(1:364))

    # TEST CHANGES TO BUDBURST ################################################################
    to_change = (DOY_Bstart = 115,)
    remSPAC_11  = remakeSPAC(discrSPAC, LAI_rel = to_change)
    # test parametrizedSPAC:
    @test remSPAC_11.parametrizedSPAC.pars.canopy_evolution.LAI_rel.DOY_Bstart != discrSPAC.parametrizedSPAC.pars.canopy_evolution.LAI_rel.DOY_Bstart
    @test remSPAC_11.parametrizedSPAC.pars.canopy_evolution.LAI_rel.DOY_Bstart == to_change.DOY_Bstart
    # test ODEProblem:
    LAI_t = remSPAC_11.ODEProblem.p.p_LAI(1:364)
    @test (to_change.DOY_Bstart) == findfirst(LAI_t .> minimum(LAI_t))

    # TEST CHANGES TO ROOTS ################################################################
    to_change = (beta = 0.88, z_rootMax_m = -0.6,)
    remSPAC_12  = remakeSPAC(discrSPAC, root_distribution = to_change)
    # test parametrizedSPAC:
    @test remSPAC_12.parametrizedSPAC.pars.root_distribution.beta        .≈ to_change.beta
    @test remSPAC_12.parametrizedSPAC.pars.root_distribution.z_rootMax_m .≈ to_change.z_rootMax_m
    # test ODEProblem:
    @test all(remSPAC_12.ODEProblem.p.p_fT_RELDEN.itp.coefs[1,:] .≈
                LWFBrook90.Rootden_beta_(
                    to_change.beta,
                    Δz_m = remSPAC_12.parametrizedSPAC.soil_discretization.Δz,
                    z_rootMax_m = to_change.z_rootMax_m))

    # code to easily modify:
    # - θs, Ks, α                    # p_THSAT, p_KSAT, p_MvGα,              all in parametrizedSPAC.pars.soil_horizons.shp[1] and proportionally all other layers
    # - drain, bypar, ilayer, infexp # DRAIN, BYPAR, BYPAR, IDEPTH_m, INFEXP all in parametrizedSPAC.pars.params
    # - glmax, r5, cvpd              # GLMAX, R5, CVPD                       all in parametrizedSPAC.pars.params
    # - cintrl                       # CINTRL                                all in parametrizedSPAC.pars.params
    # - alb, albsn                   # ALB, ALBSN                            all in parametrizedSPAC.pars.params
    # - rssa                         # RSSA                                  all in parametrizedSPAC.pars.params
    # - budburstdoy (i.e. LAI(t))    # DOY_Bstart, DOY_Cstart, LAI_perc_CtoB all in parametrizedSPAC.pars.canopy_evolution.LAI_rel
    # - psicr, fxylem, mxkpl         # PSICR, FXYLEM, MXKPL                  all in parametrizedSPAC.pars.params
    # - maxlai                       # MAXLAI                                all in parametrizedSPAC.pars.params
    # - maxrootdepth, betaroot       # beta, z_rootMax_m                     all in parametrizedSPAC.pars.root_distribution
end

@testset "simulate-and-postprocess" begin
    Δz_m = fill(0.05, 22)
    parametrizedSPAC = loadSPAC(
            "../examples/DAV2020-full/", "DAV2020-full"; simulate_isotopes = true,
            # "../examples/DAV2020-bare-minimum/", "DAV2020-minimal"; simulate_isotopes = true,
            Δz_thickness_m = Δz_m,
            root_distribution = (beta = 0.77, z_rootMax_m = -0.5),
            IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -10.11111, delta2H_init_permil = -91.1111),
            canopy_evolution = (DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel    = 100,
                                            LAI_rel = (DOY_Bstart = 120,
                                                Bduration  = 21,
                                                DOY_Cstart = 270,
                                                Cduration  = 60,
                                                LAI_perc_BtoC = 95,
                                                LAI_perc_CtoB = 70)));

    simulation = setup(parametrizedSPAC);
    simulate!(simulation)

    # check output
    depths_to_test_mm = [100, 1000, 1200, 200, 300, 400, 150, ] # test unsorted input
    @test_logs (:warn, r"below simulation domain") LWFBrook90.get_soil_idx(simulation, depths_to_test_mm)
    idx_to_read_out = LWFBrook90.get_soil_idx(simulation, depths_to_test_mm)
    @test idx_to_read_out == Dict(
        150  => 3, 100  => 2, 200  => 4, 300  => 6, 400  => 8, 1000 => 21,
        1200 => 0) # 1200 is below the simulation domain!
    valid_idx_to_read_out = LWFBrook90.get_soil_idx(simulation, depths_to_test_mm; only_valid_idxs = true)

    depths_to_test_mm_noWarning = [100, 1000, 200, 300, 400, 150, ]
    @test_throws r"ambiguous which colum represents" all(-7 .≈ get_soil_([:ψ,:θ], simulation; depths_to_read_out_mm = depths_to_test_mm_noWarning, days_to_read_out_d = 0, flag_return_Matrix = true)) # flag_return_Matrix is deprecated
    @test all(-7            .≈ get_soil_(:ψ, simulation; depths_to_read_out_mm = depths_to_test_mm_noWarning, days_to_read_out_d = 0, flag_return_Matrix = true)) # flag_return_Matrix is deprecated
    @test all(0.20058687988 .≈ get_soil_(:θ, simulation; depths_to_read_out_mm = depths_to_test_mm_noWarning, days_to_read_out_d = 0, flag_return_Matrix = true)) # flag_return_Matrix is deprecated
    @test all(-10.11111 .≈ get_soil_(:δ18O, simulation; depths_to_read_out_mm = depths_to_test_mm_noWarning, days_to_read_out_d = 0, flag_return_Matrix = true)) # flag_return_Matrix is deprecated
    @test all(-10.11111 .≈ get_soil_(:delta18O, simulation; depths_to_read_out_mm = depths_to_test_mm_noWarning, days_to_read_out_d = 0, flag_return_Matrix = true)) # flag_return_Matrix is deprecated
    @test all(-91.1111  .≈ get_soil_(:δ2H, simulation; depths_to_read_out_mm = depths_to_test_mm_noWarning, days_to_read_out_d = 0, flag_return_Matrix = true)) # flag_return_Matrix is deprecated

    # get_soil_(:θ, simulation)
    # get_soil_([:θ], simulation)
    # get_soil_([:θ], simulation; depths_to_read_out_mm = [100, 200, 500, 1200])
    # get_soil_([:θ, :ψ], simulation; depths_to_read_out_mm = [100, 200, 500, 1200])
    # get_soil_([:θ, :ψ, :K], simulation; depths_to_read_out_mm = [100, 200, 500, 1200])
    # get_soil_([:θ, :ψ, :W, :SWATI, :K], simulation; depths_to_read_out_mm = [100, 500])
    # get_soil_([:θ, :ψ, :δ18O, :δ2H, :W, :SWATI, :K], simulation; depths_to_read_out_mm = [100, 200, 500, 1200])

    # Check water partitioning output
    @test_throws r"daily resolution" get_water_partitioning(simulation)
    simulate!(simulation; save_everystep = false, saveat = range(parametrizedSPAC.tspan...), tspan = parametrizedSPAC.tspan);
    df_partitioning_daily, df_partitioning_monthly, df_partitioning_yearly = get_water_partitioning(simulation)

    reference_daily_check = [
    0.0419232  0.0092856  0.0        0.0      0.0326376  0.0  -2.08167e-17  -0.223151  -0.0   83.8668;    # DateTime(2021-01-10)
    0.0485448  0.0        0.0166119  0.0      0.0319328  0.0   3.46945e-17  -0.217801  -0.0   83.8769;    # DateTime(2021-01-20)
    0.437416   0.0        0.215244   0.0      0.222172   0.0   0.0          -2.10676   -0.0  104.066;     # DateTime(2021-04-10)
    1.62382    0.0038898  0.0        0.2079   1.41203    0.9   1.18509      -0.071398  -0.0   55.9741;    # DateTime(2021-06-29)
    1.85095    0.180715   0.0        0.0      1.67024    0.0  -2.22045e-16  -1.0321    -0.0   91.48;      # DateTime(2021-09-07)
    1.26533    0.0277206  0.0        1.08752  0.150094   7.7  -2.77556e-17  -0.262046  -0.0   76.4632;    # DateTime(2021-09-17)
    ]
    computed_daily_check = Matrix(df_partitioning_daily[[10, 20, 100, 180, 250, 260],
        ["ETa", "Esoil", "Esnow", "Einterception", "Ta", "Precip", "Td", "D", "R", "Swat",]])
    @test isapprox(reference_daily_check, computed_daily_check, atol = 1e-5, rtol = 1e-5)

    reference_monthly_check = [
    12.6084    0.190331  1.02305    10.6697    0.725283  158.1   1.01915e-17    -6.53731   0.0   84.8226
    13.0302    0.0       3.56712     4.97      4.49312    24.5  -1.56125e-16   -44.2379    0.0  103.0
    18.3249    0.0       5.6995      6.73528   5.8901     50.4  -2.68882e-16   -50.7475    0.0  101.607
    22.7934    1.79797   3.98622     4.99301  12.0162     30.2  -9.71445e-17   -90.2158    0.0  104.712
    44.5537    5.51835   0.0        15.3928   23.6426     91.3  -8.04912e-16   -45.9208    0.0  101.613
    67.0981    2.13576   0.0        10.5567   54.4056     45.7   1.9749         -8.95216   0.0   72.0244
    65.2563    6.20792   0.0        22.5187   36.5297    159.4   5.55112e-17   -41.8129    0.0   92.1359
    58.3262    6.58114   0.0        17.0969   34.6482    181.3   4.85723e-16  -129.34      0.0  106.741
    47.1773    4.38339   0.0         8.77759  34.0164     56.2   1.19349e-15   -18.2377    0.0   88.4071
    22.4298    3.01918   0.0         2.39139  17.0192     12.1   1.09635e-15   -14.1241    0.0   86.5785
    10.8751    0.823499  0.894483    5.50742   3.64969    83.7   1.38778e-17   -50.928     0.0  100.428
    9.9458    0.0       2.21222     6.58132   1.15226    63.7  -1.26201e-16   -12.7609    0.0   89.5052
    0.170659  0.0       0.0488657   0.0       0.121793    0.0   0.0            -0.813725  0.0  113.958
    ]
    computed_monthly_check = Matrix(df_partitioning_monthly[:,
        ["ETa", "Esoil", "Esnow", "Einterception", "Ta", "Precip", "Td", "D", "R", "Swat",]])
    @test isapprox(reference_monthly_check, computed_monthly_check, atol = 1e-5, rtol = 1e-5)


    reference_yearly_check = [
        392.419     30.6575  17.3826     116.191  228.188     956.6  1.9749  -513.815     0.0   94.2582
        0.170659   0.0      0.0488657    0.0      0.121793    0.0  0.0       -0.813725  0.0  113.958
    ]
    computed_yearly_check = Matrix(df_partitioning_yearly[:,
        ["ETa", "Esoil", "Esnow", "Einterception", "Ta", "Precip", "Td", "D", "R", "Swat",]])
    @test isapprox(reference_yearly_check, computed_yearly_check, atol = 1e-5, rtol = 1e-5)
end
