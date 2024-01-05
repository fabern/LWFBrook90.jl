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

@testset "prepare SPAC for LWFBrook90" begin
    Δz = [fill(0.05, 4);  fill(0.10, 14)]
    input_path   = "../examples/DAV2020-full/";
    model        = loadSPAC(input_path, "DAV2020-full"; simulate_isotopes = false);
    mod_model    = loadSPAC(input_path, "DAV2020-full"; simulate_isotopes = false,
            Δz_thickness_m    = Δz,
            root_distribution = (beta = 0.98, z_rootMax_m = -sum(Δz)), # use whole domain as root zone (beta parameter regulates distribution)
            IC_soil           = (PSIM_init_kPa = -6.0,
                                delta18O_init_permil = -9.0,
                                delta2H_init_permil = -11.0),
            canopy_evolution  = (DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel    = 100,
                                LAI_rel = (DOY_Bstart = 110,    Bduration  = 20,
                                            DOY_Cstart = 270,    Cduration  = 60,
                                            LAI_perc_BtoC = 100, LAI_perc_CtoB = 0)),
            storm_durations_h = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
            IC_scalar         = (amount = (u_GWAT_init_mm = 0.0,         u_INTS_init_mm = 0.0,
                                        u_INTR_init_mm = 0.0,         u_SNOW_init_mm = 0.0,
                                        u_CC_init_MJ_per_m2 = 0.0001, u_SNOWLQ_init_mm =  0.),
                                d18O    = (u_GWAT_init_permil = -11.111, u_INTS_init_permil = -12.222,
                                        u_INTR_init_permil = -13.333, u_SNOW_init_permil = -14.444),
                                d2H     = (u_GWAT_init_permil = -95.111, u_INTS_init_permil = -95.222,
                                        u_INTR_init_permil = -95.333, u_SNOW_init_permil = -95.444)));

    base_simulation = LWFBrook90.setup(model; requested_tspan = (0,300));
    mod_simulation  = LWFBrook90.setup(mod_model)


    # NOTE: in a next iteration, we could support LAI as vector (among other stuff) by
    #   generating input data for r_lwfbrook90() instead of (effectively skipping preprocessing
    #   with the R-pkg and directly going for the Fortran Code.)
    @test_throws r"LAI as vector is not supported"      args, derived_args = prepare_for_LWFBrook90R(base_simulation, return_value = "inputs");
    @test_throws r"Please provide a `discretized` SPAC" args, derived_args = prepare_for_LWFBrook90R(mod_model, return_value = "inputs");

    args, derived_args = prepare_for_LWFBrook90R(mod_simulation, return_value = "inputs");

    # b1)
    @test isequal(Matrix(args.meteo[365:366,:]),
        #generated with: print(IOContext(stdout, :compact=>false), Matrix(args.meteo[365:366,:]))
        Any[Date("2021-12-31") -1.8 6.6 NaN 0.0 100 6.13 1.1 0.52; Date("2022-01-01") -1.8 6.6 NaN 0.0 100 6.13 1.1 0.52])
    @test isapprox(
        Matrix(args.soil[:, [:upper, :lower, :gravel, :ths, :thr, :alpha, :npar, :mpar, :ksat, :tort, :rootden]]),
        [       # generated with print(IOContext(stdout, :compact=>false), Matrix(args.soil[:, [:upper, :lower, :gravel, :ths, :thr, :alpha, :npar, :mpar, :ksat, :tort, :rootden]]))
        -0.0 -0.05 0.375 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.14862236000551818;
        -0.05 -0.1 0.375 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.1343428420784842;
        -0.1 -0.15000000000000002 0.375 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.12143528885596006;
        -0.15000000000000002 -0.2 0.375 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.10976788306231737;
        -0.2 -0.3 0.375 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.09445491232028565;
        -0.3 -0.39999999999999997 0.375 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.07717654033385292;
        -0.39999999999999997 -0.49999999999999994 0.75 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.06305885243645104;
        -0.49999999999999994 -0.6 0.75 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.05152367355935886;
        -0.6 -0.7 0.75 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.04209859257630285;
        -0.7 -0.7999999999999999 0.75 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.03439761520233502;
        -0.7999999999999999 -0.8999999999999999 0.75 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.028105356003609673;
        -0.8999999999999999 -0.9999999999999999 0.75 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.02296412211844314;
        -0.9999999999999999 -1.0999999999999999 0.75 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.018763359717024704;
        -1.0999999999999999 -1.2 0.75 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.015331030990630084;
        -1.2 -1.3 0.75 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.012526568523994051;
        -1.3 -1.4000000000000001 0.75 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.010235118504569107;
        -1.4000000000000001 -1.5000000000000002 0.75 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.008362837005354907;
        -1.5000000000000002 -1.6000000000000003 0.75 0.3786 0.0 20.387 1.2347 0.19008666072730207 2854.91 -3.339 0.006833046705508307])
    @test isequal(args.opts,
        (startdate = DateTime("2021-01-01T00:00:00"), enddate = DateTime("2021-12-31T00:00:00"),
        fornetrad = "globrad", prec_interval = 1, correct_prec = false, budburst_method = "fixed",
        leaffall_method = "fixed", standprop_input = "parameters", standprop_interp = "constant",
        use_growthperiod = false, lai_method = "b90", imodel = "MvG", root_method = "soilvar"))
    @test isequal(args.parms,
        (maxlai = 3.0, sai = 1.0, sai_ini = 1.0, height = 25.0, height_ini = 25.0, densef = 1.0,
        densef_ini = 1.0, age_ini = 100.0, winlaifrac = 0.0, budburst_species = missing,
        budburstdoy = 110, leaffalldoy = 270, emergedur = 20, leaffalldur = 60,
        shp_optdoy = missing, shp_budburst = missing, shp_leaffall = missing, alb = 0.2,
        albsn = 0.5, ksnvp = 0.3, fxylem = 0.5, mxkpl = 15.64, lwidth = 0.1, psicr = -1.03942,
        nooutf = 1, lpc = 4.0, cs = 0.35, czs = 0.13, czr = 0.05, hs = 1.0, hr = 10.0,
        rhotp = 2.0, nn = 2.5, maxrlen = 3000.0, initrlen = 12.0, initrdep = 0.25,
        rrad = 0.35, rgrorate = 0.03, rgroper = 30.0, maxrootdepth = missing, betaroot = missing,
        radex = 0.5, glmax = 0.00868, glmin = 0.0003, rm = 1000.0, r5 = 287.0,
        cvpd = 2.0, tl = 0.0, t1 = 10.0, t2 = 30.0, th = 40.0, frintlai = 0.06, frintsai = 0.06,
        fsintlai = 0.04, fsintsai = 0.04, cintrl = 0.15, cintrs = 0.15, cintsl = 0.6, cintss = 0.6,
        infexp = 0.45, bypar = 0, qfpar = 1.0, qffc = 0.0, imperv = 0.0,
        drain = 0.51, gsc = 0.0, gsp = 0.0, ilayer = 6, qlayer = 0, z0s = 0.001, rstemp = -0.5,
        melfac = 1.5, ccfac = 0.3, laimlt = 0.2, saimlt = 0.5, grdmlt = 0.35, maxlqf = 0.05,
        snoden = 0.3, obsheight = 0.024999999999999998, correct_prec_statexp = "mg",
        rssa = 795.29579, rssb = 1.0, dtimax = 0.5, dswmax = 0.05, dpsimax = 0.0005,
        wndrat = 0.3, fetch = 5000.0, z0w = 0.005, zw = 2.0, zminh = 2.0, coords_x = 0,
        coords_y = 47.0, c1 = 0.25, c2 = 0.5, c3 = 0.2,
        pdur = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
        eslope = 0.0, aspect = 0.0, dslope = 0.0,
        slopelen = 200.0, intrainini = 0.0, intsnowini = 0.0, gwatini = 0.0, snowini = 0.0,
        psiini = -6.0, standprop_table = missing, lai_doy = missing, lai_frac = missing,
        rootden_table = missing, soil_nodes = missing, soil_materials = missing))
    # # b2)
    # derived_args.meteo
    # derived_args.soil
    # derived_args.df_opts
    # derived_args.df_parms_a
    # derived_args.df_parms_b
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

    # Check isotope output for plotting
    df_isotopePlot, RWUPlotlabel = LWFBrook90.get_data_for_isotopePlot(simulation)
    @test RWUPlotlabel == "mean RWU depth\n(based on uptake only)"
    # print(IOContext(stdout, :compact=>false), permutedims(Matrix(df_isotopePlot[[1, 10, 100, 300], :])))
    reference_isotopePlotPermutedDims = [
    0.0 9.0 99.0 299.0;
    0.0 0.0 0.0 0.0;
    -15.04 -15.04 -14.42 -22.6;
    -111.96 -111.96 -102.2 -171.53;
    -15.04 -15.04 -14.42 -22.6;
    -111.96 -111.96 -102.2 -171.53;
    -13.0 -15.04 -14.42 -22.6;
    -13.0 -15.04 -14.42 -22.6;
    -13.0 -15.04 -14.42 -22.6;
    -13.0 -10.609721027648131 -14.050570887472329 -9.191276335721108;
    -10.11111 -10.261397637901904 -14.345732332233464 -10.410159261495577;
    -10.11111 -10.111715283165307 -12.271846832576891 -10.205996912280122;
    -95.0 -111.96 -102.2 -171.53;
    -95.0 -111.96 -102.2 -171.53;
    -95.0 -111.96 -102.2 -171.53;
    -95.0 -91.78230591535274 -100.66664640592148 -60.91149805549817;
    -91.1111 -91.746331679754 -101.64448413109162 -71.98553453622014;
    -91.1111 -91.11365913620521 -97.4679731838423 -69.78144021861249;
    NaN 43.22335331429749 43.2198758960631 50.39351555370525;
    -10.11111 -10.288628249440086 -14.35367998200292 -10.453980436189958;
    -10.11111 -10.195345054215933 -14.332944610392223 -10.368919551316301;
    -10.11111 -10.166563235529507 -14.309103005246701 -10.273116212689624;
    -10.11111 -10.155331406479542 -14.283926629937943 -10.190657282276014;
    -10.11111 -10.14892149017908 -14.258268733634052 -10.117169655501959;
    -10.11111 -10.144497776974289 -14.233060989383596 -10.04798302714493;
    -10.11111 -10.140913662489963 -14.209676871001333 -9.980880216128392;
    -10.11111 -10.135715781921455 -14.190260645002459 -9.913844192660417;
    -10.11111 -10.122751529843436 -14.17716820911132 -9.843925427630435;
    -10.11111 -10.11491223862554 -14.164270950588909 -9.778612380690804;
    -10.11111 -10.112052211355714 -14.151405237844935 -9.714618777485288;
    -10.11111 -10.111298004185059 -14.138627395482938 -9.651317903636167;
    -10.11111 -10.111141360898204 -14.126098081307303 -9.589483325769578;
    -10.11111 -10.111114482599989 -14.114001479130213 -9.530076065477164;
    -10.11111 -10.111110558420712 -14.102514057214119 -9.47323198479031;
    -10.11111 -10.1111100613681 -14.091800486771957 -9.419140109014062;
    -10.11111 -10.11111000600247 -14.08201596012743 -9.36853253672173;
    -10.11111 -10.111110000526027 -14.073304614300728 -9.32232462973745;
    -10.11111 -10.111110000041515 -14.065795008604228 -9.281618649114554;
    -10.11111 -10.111110000002947 -14.059606912856875 -9.247803312939103;
    -10.11111 -10.111110000000174 -14.054904158286385 -9.222700581645364;
    -10.11111 -10.11111000000004 -14.05206217265466 -9.208688356421478;
    -91.1111 -91.86117687621623 -101.70088056934627 -72.41648152042562;
    -91.1111 -91.46794441874937 -101.55187485661764 -71.58588734686599;
    -91.1111 -91.34582429241216 -101.38589001517572 -70.63189647388069;
    -91.1111 -91.29819728949543 -101.21707922140061 -69.80244614854549;
    -91.1111 -91.27105694261236 -101.05222655204668 -69.06971706146682;
    -91.1111 -91.25233453564253 -100.89749731844239 -68.39238337588111;
    -91.1111 -91.23715180777555 -100.76039710907169 -67.74797747236411;
    -91.1111 -91.21510565408296 -100.65155242466521 -67.115809209937;
    -91.1111 -91.16046275330156 -100.58189524654347 -66.46806297726769;
    -91.1111 -91.12729530639031 -100.51830882239275 -65.86721978050565;
    -91.1111 -91.11513373815411 -100.4606710756425 -65.28404977381996;
    -91.1111 -91.11190921633201 -100.41042734587118 -64.71387158903909;
    -91.1111 -91.11123573854987 -100.36980306837496 -64.1647849350807;
    -91.1111 -91.1111195128364 -100.34122158039206 -63.646080890953876;
    -91.1111 -91.11110244497971 -100.32702267562412 -63.159084129549264;
    -91.1111 -91.1111002702853 -100.32929790500059 -62.705081808694366;
    -91.1111 -91.11110002659602 -100.34964449568014 -62.28940535002614;
    -91.1111 -91.1111000023449 -100.38865988340882 -61.91810943103562;
    -91.1111 -91.11110000018614 -100.44498717419589 -61.597892041300796;
    -91.1111 -91.11110000001332 -100.51368501758468 -61.336906141380624;
    -91.1111 -91.11110000000085 -100.58372491875083 -61.14611169397025;
    -91.1111 -91.1111 -100.6347904184193 -61.04068884469542]

    @testset "postprocess_isotope-plot" begin
        [@test filter(x -> !isnan(x), c) ≈ filter(x -> !isnan(x), r)
            for ((k1, c), (k2, r)) in
            zip(pairs(eachcol(df_isotopePlot[[1, 10, 100, 300], :])),
                pairs(eachcol(permutedims(reference_isotopePlotPermutedDims))))]
    end

end
