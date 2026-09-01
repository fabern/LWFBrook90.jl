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

@testset "Parameter defaults and incomplete files" begin
    defaults = LWFBrook90.default_param_values()
    @test all(isnan, [defaults.LAT_DEG, defaults.ESLOPE_DEG, defaults.ASPECT_DEG, defaults.MAXLAI])
    @test defaults.DTIMAX == 0.5

    # test use of incomplete param.csv file, where missing parameters can be filled with defaults
    # test use of incomplete param.csv file, where missing parameters can NOT be filled with defaults
    path, io = mktemp()
    try
        Base.write(io, "param_id,x\nLAT_DEG,48.5\nESLOPE_DEG,10\nASPECT_DEG,180\nMAXLAI,4\nNOOUTF,0\n")
        close(io)

        params, solver_options = LWFBrook90.read_path_param(path)
        @test params.LAT_DEG == 48.5
        @test params.ESLOPE_DEG == 10.0
        @test params.ASPECT_DEG == 180.0
        @test params.MAXLAI == 4.0
        @test params.NOOUTF == 0
        @test params.ALB == defaults.ALB
        # and others in between ...
        @test solver_options.DTIMAX == defaults.DTIMAX
        @test solver_options.DSWMAX == defaults.DSWMAX
        @test solver_options.DPSIMAX == defaults.DPSIMAX
    finally
        isopen(io) && close(io)
        rm(path; force = true)
    end

    # test use of incomplete param.csv file, where missing parameters can NOT be filled with defaults
    path, io = mktemp()
    try
        Base.write(io, "param_id,x\nLAT_DEG,48.5\nESLOPE_DEG,NaN\n")
        close(io)

        thrown_error = try
            LWFBrook90.read_path_param(path)
            nothing
        catch error
            error
        end
        @test thrown_error isa ErrorException
        @test occursin("Missing or NaN: [\"ESLOPE_DEG\", \"ASPECT_DEG\", \"MAXLAI\"]", thrown_error.msg)
    finally
        isopen(io) && close(io)
        rm(path; force = true)
    end
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
    f1 = (Δz_m) -> LWFBrook90.Rootden_(beta = 0.97, Δz_m = Δz_m)  # function for root density as f(Δz)
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

@testset "root-model gamma (Δz_m = $(first(Δz_m)))" for Δz_m in Δz_m_data # source: https://stackoverflow.com/a/63871951
    # Check equivalence of beta model and gamma model (with k=1.0)
    @test LWFBrook90.Rootden_(; beta = 0.97, Δz_m = Δz_m) ≈
        LWFBrook90.Rootden_(; root_k = 1.0, root_θ_cm = -1/log(0.97), Δz_m = Δz_m)
    @test LWFBrook90.Rootden_(; beta = 0.90, Δz_m = Δz_m) ≈
        LWFBrook90.Rootden_(; root_k = 1.0, root_θ_cm = -1/log(0.90), Δz_m = Δz_m)

    # check
    parametrizedSPAC = loadSPAC("test-assets/Hammel-2001/input-files-ISO",
                      "Hammel_loam-NLayer-27-RESET=FALSE";
                      simulate_isotopes = false,
                      Δz_thickness_m = Δz_m,
                      root_distribution = (root_k = 1.0, root_θ_cm = -1/log(0.95),),
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

    # checks with k>1.0
    parametrizedSPAC = loadSPAC("test-assets/Hammel-2001/input-files-ISO",
                      "Hammel_loam-NLayer-27-RESET=FALSE";
                      simulate_isotopes = false,
                      Δz_thickness_m = Δz_m,
                      root_distribution = (root_k = 3.0, root_θ_cm = -1/log(0.95),),
                      IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -9.0, delta2H_init_permil = -11.0))
    # parametrizedSPAC.pars.root_distribution
    simulation = setup(parametrizedSPAC, ε = 0.005)
    # plot(simulation.parametrizedSPAC.soil_discretization.df.Rootden_, simulation.parametrizedSPAC.soil_discretization.df.Lower_m)
    if ([0.5, 0.5, 0.5, 0.5] == Δz_m)
        # @test simulation.parametrizedSPAC.soil_discretization.df.Rootden_ ≈ [0.018461100494465737, 0.018461100494465737, 0.0014204889211275828, 0.00010929948491700703, 8.410046164697427e-6]
        @test simulation.parametrizedSPAC.soil_discretization.df.Rootden_ ≈ [0.00027502637688857145, 0.47609855754620745, 0.412005500758867, 0.09645815471805702, 0.015162760599979806]
    elseif ([fill(0.02, 100);] == Δz_m)
        @test simulation.parametrizedSPAC.soil_discretization.df.Rootden_[1:10:100] ≈ [1.105600935589554e-5, 0.01841592534143428, 0.027821629343506862, 0.02282995046421708, 0.014674812772489988, 0.008262099744807986, 0.0042796350631698135, 0.0020932916762929975, 0.0009819246729527453, 0.0004461401605346446]
    elseif ([fill(0.01, 200);] == Δz_m)
        @test simulation.parametrizedSPAC.soil_discretization.df.Rootden_[1:10:200] ≈ [5.528233262478093e-6, 0.003749688650354093, 0.009456878548333922, 0.012958944185572032, 0.013911350780858728, 0.013080734304564738, 0.011316177973638309, 0.009244379118960236, 0.007242414922077211, 0.005495855155572414, 0.004067005209129574, 0.0029491416245713367, 0.002103011304925199, 0.0014787113599272316, 0.0010273776115250128, 0.0007064829092284812, 0.0004814797904216458, 0.0003255617365651601, 0.00021860495106492463, 0.00014587690668812743]
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
    
    @test_throws r"First discretize it with .*setup\(SPAC\)" simulate!(parametrizedSPAC)
    
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
    @test all(remSPAC_12.ODEProblem.p.p_RELDEN.itp.coefs[1,:] .≈
                LWFBrook90.Rootden_(
                    beta = to_change.beta,
                    Δz_m = remSPAC_12.parametrizedSPAC.soil_discretization.Δz,
                    z_rootMax_m = to_change.z_rootMax_m))

    to_change = (root_θ_cm = 20, root_k = 3.0, z_rootMax_m = -0.93)
    remSPAC_12  = remakeSPAC(discrSPAC, root_distribution = to_change)

    # test parametrizedSPAC:
    @test remSPAC_12.parametrizedSPAC.pars.root_distribution.root_θ_cm   .≈ to_change.root_θ_cm
    @test remSPAC_12.parametrizedSPAC.pars.root_distribution.root_k      .≈ to_change.root_k
    @test remSPAC_12.parametrizedSPAC.pars.root_distribution.z_rootMax_m .≈ to_change.z_rootMax_m
    # test ODEProblem:
    must = LWFBrook90.Rootden_(
                    root_k = to_change.root_k, root_θ_cm = to_change.root_θ_cm,
                    Δz_m = remSPAC_12.parametrizedSPAC.soil_discretization.Δz,
                    z_rootMax_m = to_change.z_rootMax_m)
    must = [0.0019726404225474856, 0.023343292699885786, 0.07757036883704736, 0.13045558629384227, 0.15554436815528225, 0.1563597889663928, 0.1419548680219542, 0.12044192366907798, 0.09737691045459235, 0.07593633493301175, 0.019043917546365884, 0.0]
    findmax(must)[2] == 6
    @test all(remSPAC_12.ODEProblem.p.p_RELDEN.itp.coefs[1,:] .≈ must)

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
        Any[Date("2021-12-31") -1.8 6.6 NaN 0.0 100 6.13 1.1 0.52; 
            Date("2022-01-01") -1.8 6.6 NaN 0.0 100 6.13 1.1 0.52])
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

@testset "check consistency of output (API)" begin
    # case 1 simulate_isotopes = FALSE, simulate_irrigation = FALSE
    # case 2 simulate_isotopes = TRUE, simulate_irrigation = FALSE
    # case 3 simulate_isotopes = FALSE, simulate_irrigation = TRUE
    # case 4 simulate_isotopes = TRUE, simulate_irrigation = TRUE

    # setup and run simulation
    parametrizedSPAC = loadSPAC(
            "../examples/DAV2020-full/", "DAV2020-full"; 
            simulate_isotopes = true, simulate_irrigation = false,
            Δz_thickness_m = fill(0.05, 22),
            root_distribution = (beta = 0.77, z_rootMax_m = -0.5),
            IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -10.11111, delta2H_init_permil = -91.1111),
            canopy_evolution = (
                DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel = 100,
                LAI_rel = (DOY_Bstart = 120,
                           Bduration  = 21,
                           DOY_Cstart = 270,
                           Cduration  = 60,
                           LAI_perc_BtoC = 95,
                           LAI_perc_CtoB = 70)));

    simulation = setup(parametrizedSPAC);
    simulate!(simulation)

    # Check requirement of daily output resolution for some outputs:
    @test_throws r"Solution is not computed with daily output resolution." get_water_partitioning(simulation)
    # TODO: do we need that check for get_states, get_fluxes as well??

    # Rerun with daily:
    simulate!(simulation; 
        save_everystep = false, 
        saveat = range(parametrizedSPAC.tspan...), 
        tspan = parametrizedSPAC.tspan);

    # A) check consistency of documented output functions: 
    # (get_states, get_fluxes, get_forcing, get_soil_, get_water_partitioning)
    depths_to_test_mm_noWarning = [100, 1000, 200, 300, 400, 150, ]

    # compute current output to test
    df_get_forcing = get_forcing(simulation)[[1],:] # must be a DataFrame not a DataFrameRow
    df_get_states  = get_states(simulation, days_to_read_out_d = 0)
    df_get_fluxes  = get_fluxes(simulation, days_to_read_out_d = 0)
    df_get_soil_all= get_soil_([:θ,:ψ,:W,:SWATI,:K,:δ18O,:TRANI],
                               simulation, depths_to_read_out_mm = depths_to_test_mm_noWarning, days_to_read_out_d = 0)
    df_partitioning_daily, df_partitioning_monthly, df_partitioning_yearly = 
        get_water_partitioning(simulation)

    # get reference values of DataFrames (not DataFrameRows !) with: 
    function julian_dput(df::DataFrame; multiline = false)
        sep = ifelse(multiline, "\n","")
        cols = [repr(df[!, col]) for col in names(df)]
        code = "DataFrame($sep$(join(["$c=$v" for (c,v) in zip(names(df), cols)], ",$sep")))"
        println(code)
        return nothing
    end

    # julian_dput(df_get_forcing, multiline=true)
    ref_get_forcing  = DataFrame(
        dates=[DateTime("2021-01-01T00:00:00")],
        globrad_MJDayM2=[5.53],
        tmax_degC=[0.9],
        tmin_degC=[-10.1],
        vappres_kPa=[0.26],
        windspeed_ms=[1.5],
        prec_mmDay=[0.0],
        precdelta18O_permil=[-15.04],
        precdelta2H_permil=[-111.96],
        irrig_mmDay=[0.0],
        irrigdelta18O_permil=[nothing],
        irrigdelta2H_permil=[nothing],
        densef_percent=[100.0],
        height_m=[25.0],
        lai_m2m2=[2.0999999999999996],
        sai_m2m2=[1.0])

    # julian_dput(df_get_states, multiline=true)
    ref_get_states   = DataFrame(
        dates=[DateTime("2021-01-01T00:00:00")],
        GWAT_mm=[1.0],
        INTS_mm=[0.0],
        INTR_mm=[0.0],
        SNOW_mm=[0.0],
        CC_MJm2=[0.0],
        SNOWLQ_mm=[0.0],
        SWAT_mm=[85.24942394707108],
        GWAT_d18O=[-13.0],
        INTS_d18O=[-13.0],
        INTR_d18O=[-13.0],
        SNOW_d18O=[-13.0],
        XYLEM_d18O=[-10.11111],
        GWAT_d2H=[-95.0],
        INTS_d2H=[-95.0],
        INTR_d2H=[-95.0],
        SNOW_d2H=[-95.0],
        XYLEM_d2H=[-91.1111],
        d18O_permil_50mm=[-10.11111],
        d18O_permil_100mm=[-10.11111],
        d18O_permil_150mm=[-10.11111],
        d18O_permil_200mm=[-10.11111],
        d18O_permil_250mm=[-10.11111],
        d18O_permil_300mm=[-10.11111],
        d18O_permil_350mm=[-10.11111],
        d18O_permil_400mm=[-10.11111],
        d18O_permil_450mm=[-10.11111],
        d18O_permil_500mm=[-10.11111],
        d18O_permil_550mm=[-10.11111],
        d18O_permil_600mm=[-10.11111],
        d18O_permil_650mm=[-10.11111],
        d18O_permil_700mm=[-10.11111],
        d18O_permil_750mm=[-10.11111],
        d18O_permil_800mm=[-10.11111],
        d18O_permil_850mm=[-10.11111],
        d18O_permil_900mm=[-10.11111],
        d18O_permil_950mm=[-10.11111],
        d18O_permil_1000mm=[-10.11111],
        d18O_permil_1050mm=[-10.11111],
        d18O_permil_1100mm=[-10.11111],
        d2H_permil_50mm=[-91.1111],
        d2H_permil_100mm=[-91.1111],
        d2H_permil_150mm=[-91.1111],
        d2H_permil_200mm=[-91.1111],
        d2H_permil_250mm=[-91.1111],
        d2H_permil_300mm=[-91.1111],
        d2H_permil_350mm=[-91.1111],
        d2H_permil_400mm=[-91.1111],
        d2H_permil_450mm=[-91.1111],
        d2H_permil_500mm=[-91.1111],
        d2H_permil_550mm=[-91.1111],
        d2H_permil_600mm=[-91.1111],
        d2H_permil_650mm=[-91.1111],
        d2H_permil_700mm=[-91.1111],
        d2H_permil_750mm=[-91.1111],
        d2H_permil_800mm=[-91.1111],
        d2H_permil_850mm=[-91.1111],
        d2H_permil_900mm=[-91.1111],
        d2H_permil_950mm=[-91.1111],
        d2H_permil_1000mm=[-91.1111],
        d2H_permil_1050mm=[-91.1111],
        d2H_permil_1100mm=[-91.1111],
        θ_m3m3_50mm=[0.20058687987546137],
        θ_m3m3_100mm=[0.20058687987546137],
        θ_m3m3_150mm=[0.20058687987546137],
        θ_m3m3_200mm=[0.20058687987546137],
        θ_m3m3_250mm=[0.20058687987546137],
        θ_m3m3_300mm=[0.20058687987546137],
        θ_m3m3_350mm=[0.20058687987546137],
        θ_m3m3_400mm=[0.20058687987546137],
        θ_m3m3_450mm=[0.20058687987546137],
        θ_m3m3_500mm=[0.20058687987546137],
        θ_m3m3_550mm=[0.20058687987546137],
        θ_m3m3_600mm=[0.20058687987546137],
        θ_m3m3_650mm=[0.20058687987546137],
        θ_m3m3_700mm=[0.20058687987546137],
        θ_m3m3_750mm=[0.20058687987546137],
        θ_m3m3_800mm=[0.20058687987546137],
        θ_m3m3_850mm=[0.20058687987546137],
        θ_m3m3_900mm=[0.20058687987546137],
        θ_m3m3_950mm=[0.20058687987546137],
        θ_m3m3_1000mm=[0.2005868798754614],
        θ_m3m3_1050mm=[0.20058687987546137],
        θ_m3m3_1100mm=[0.20058687987546137],
        ψ_kPa_50mm=[-7.0000000000000036],
        ψ_kPa_100mm=[-7.0000000000000036],
        ψ_kPa_150mm=[-7.0000000000000036],
        ψ_kPa_200mm=[-7.0000000000000036],
        ψ_kPa_250mm=[-7.0000000000000036],
        ψ_kPa_300mm=[-7.0000000000000036],
        ψ_kPa_350mm=[-7.0000000000000036],
        ψ_kPa_400mm=[-7.0000000000000036],
        ψ_kPa_450mm=[-7.0000000000000036],
        ψ_kPa_500mm=[-7.0000000000000036],
        ψ_kPa_550mm=[-7.0000000000000036],
        ψ_kPa_600mm=[-7.0000000000000036],
        ψ_kPa_650mm=[-7.0000000000000036],
        ψ_kPa_700mm=[-7.0000000000000036],
        ψ_kPa_750mm=[-7.0000000000000036],
        ψ_kPa_800mm=[-7.0000000000000036],
        ψ_kPa_850mm=[-7.0000000000000036],
        ψ_kPa_900mm=[-7.0000000000000036],
        ψ_kPa_950mm=[-7.0000000000000036],
        ψ_kPa_1000mm=[-6.9999999999999964],
        ψ_kPa_1050mm=[-7.0000000000000036],
        ψ_kPa_1100mm=[-7.0000000000000036])
    
    # julian_dput(df_get_fluxes, multiline=true)
    ref_get_fluxes   = DataFrame(
        dates=[DateTime("2021-01-01T00:00:00")],
        cum_d_prec=[0.0],
        cum_d_sfal=[0.0],
        cum_d_sthr=[0.0],
        cum_d_sint=[0.0],
        cum_d_irrig=[0.0],
        cum_d_rfal=[0.0],
        cum_d_rint=[0.0],
        cum_d_rthr=[0.0],
        cum_d_rsno=[0.0],
        cum_d_rnet=[0.0],
        cum_d_smlt=[0.0],
        cum_d_irvp=[0.0],
        cum_d_isvp=[0.0],
        cum_d_snvp=[0.0],
        cum_d_slvp=[0.04482920026138777],
        cum_d_tran=[0.056616161386777504],
        RWU=[0.056616161386777504],
        INTS=[0.0],
        INTR=[0.0],
        SNOW=[0.0],
        CC=[0.0],
        SNOWLQ=[0.0],
        RWU_d18O=[-10.11111],
        RWU_d2H=[-91.1111],
        PREC_d18O=[-15.04],
        PREC_d2H=[-111.96],
        cum_d_pint=[2.6465243861228247],
        cum_d_ptran=[0.056616161386777504],
        cum_d_pslvp=[0.36472665495887147],
        srfl=[0.0],
        slfl=[0.0],
        byfl=[0.0],
        dsfl=[0.0],
        gwfl=[0.15412189586625977],
        vrfln=[0.15412189586625977],
        flow=[0.15412189586625977],
        seep=[0.0],
        evap=[0.10144536164816528],
        TRANI_mmday_50mm=[0.043498767636911864],
        TRANI_mmday_100mm=[0.009959724394658772],
        TRANI_mmday_150mm=[0.002409761574167726],
        TRANI_mmday_200mm=[0.0005748250905928008],
        TRANI_mmday_250mm=[0.0001346303525887302],
        TRANI_mmday_300mm=[3.0767571641898455e-5],
        TRANI_mmday_350mm=[6.792328768255243e-6],
        TRANI_mmday_400mm=[1.4228398921497332e-6],
        TRANI_mmday_450mm=[-5.304024446942134e-7],
        TRANI_mmday_500mm=[0.0],
        TRANI_mmday_550mm=[0.0],
        TRANI_mmday_600mm=[0.0],
        TRANI_mmday_650mm=[0.0],
        TRANI_mmday_700mm=[0.0],
        TRANI_mmday_750mm=[0.0],
        TRANI_mmday_800mm=[0.0],
        TRANI_mmday_850mm=[0.0],
        TRANI_mmday_900mm=[0.0],
        TRANI_mmday_950mm=[0.0],
        TRANI_mmday_1000mm=[0.0],
        TRANI_mmday_1050mm=[0.0],
        TRANI_mmday_1100mm=[0.0])
    
    # julian_dput(df_get_soil_all, multiline=true)
    ref_get_soil_all = DataFrame(
        time=[0], # TODO: note this is time as Int
        K_mmday_100mm=[1.1082021409141438],
        K_mmday_150mm=[1.1082021409141438],
        K_mmday_200mm=[1.1082021409141438],
        K_mmday_300mm=[1.1082021409141438],
        K_mmday_400mm=[1.1082021409141438],
        K_mmday_1000mm=[1.1082021409141432],
        SWATI_mm_100mm=[6.268339996108168],
        SWATI_mm_150mm=[6.268339996108169],
        SWATI_mm_200mm=[6.268339996108166],
        SWATI_mm_300mm=[6.268339996108166],
        SWATI_mm_400mm=[6.268339996108166],
        SWATI_mm_1000mm=[2.507335998443264],
        TRANI_mmday_100mm=[0.0],
        TRANI_mmday_150mm=[0.0],
        TRANI_mmday_200mm=[0.0],
        TRANI_mmday_300mm=[0.0],
        TRANI_mmday_400mm=[0.0],
        TRANI_mmday_1000mm=[0.0],
        W___100mm=[0.5298121496974679],
        W___150mm=[0.5298121496974679],
        W___200mm=[0.5298121496974679],
        W___300mm=[0.5298121496974679],
        W___400mm=[0.5298121496974679],
        W___1000mm=[0.529812149697468],
        δ18O_permil_100mm=[-10.11111],
        δ18O_permil_150mm=[-10.11111],
        δ18O_permil_200mm=[-10.11111],
        δ18O_permil_300mm=[-10.11111],
        δ18O_permil_400mm=[-10.11111],
        δ18O_permil_1000mm=[-10.11111],
        θ_m3m3_100mm=[0.20058687987546137],
        θ_m3m3_150mm=[0.20058687987546137],
        θ_m3m3_200mm=[0.20058687987546137],
        θ_m3m3_300mm=[0.20058687987546137],
        θ_m3m3_400mm=[0.20058687987546137],
        θ_m3m3_1000mm=[0.2005868798754614],
        ψ_kPa_100mm=[-7.0000000000000036],
        ψ_kPa_150mm=[-7.0000000000000036],
        ψ_kPa_200mm=[-7.0000000000000036],
        ψ_kPa_300mm=[-7.0000000000000036],
        ψ_kPa_400mm=[-7.0000000000000036],
        ψ_kPa_1000mm=[-6.9999999999999964])

    # compare current and reference output
    # set up approximate comparison of DataFrames (including nothgin, DateTime, ...)
    # https://discourse.julialang.org/t/is-there-a-package-to-compare-if-two-dataframes-are-the-same/44381/4
    mycompare(x::Nothing, y::Nothing; kwargs...) = true
    mycompare(x::Number, y::Number; rtol::Real=√eps(), atol::Real=0) = isapprox(x, y)
    mycompare(x, y; kwargs...) = isequal(x, y)
    function myDFcompare(x::DataFrame, y::DataFrame;
                      rtol::Real=√eps(), atol::Real=0) 
        print(x);print("\n")
        print(y);print("\n")
        all([all(vec) for vec in eachcol(mycompare.(x, y))])
    end

    @test myDFcompare(df_get_forcing, ref_get_forcing)
    @test myDFcompare(df_get_forcing, ref_get_forcing)
    @test myDFcompare(df_get_states,  ref_get_states)
    @test myDFcompare(df_get_fluxes,  ref_get_fluxes)
    # @test myDFcompare(df_get_soil_1,  ref_get_soil_1)
    # @test myDFcompare(df_get_soil_2,  ref_get_soil_2)
    # @test myDFcompare(df_get_soil_3,  ref_get_soil_3)
    # @test myDFcompare(df_get_soil_4,  ref_get_soil_4)
    # @test myDFcompare(df_get_soil_5,  ref_get_soil_5)
    # @test myDFcompare(df_get_soil_6,  ref_get_soil_6)
    # @test myDFcompare(df_get_soil_7,  ref_get_soil_7)

    # Check output of water partitioning output
    cols_to_check = ["ETa", "Esoil", "Esnow", "Einterception", "Ta", "Precip", "Td", "D", "R", "Swat",]
    # julian_dput(df_partitioning_daily[[10, 20, 100, 180, 250, 260], cols_to_check], multiline=true)
    ref_partitioning_daily = DataFrame(
        ETa           = [0.041923240612276605, 0.04854475293069882, 0.4374160550376799, 1.6240215401844946, 1.8509559669497675, 1.2653333329108203],
        Esoil         = [0.009285599165961049, 0.0, 0.0, 0.0038900536174205622, 0.18071571762165603, 0.027720901826861828],
        Esnow         = [0.0, 0.016611946626641515, 0.21524394369370514, 0.0, 0.0, 0.0],
        Einterception = [0.0, 0.0, 0.0, 0.2079, 0.0, 1.0875185605338626],
        Ta            = [0.032637641446315556, 0.0319328063040573, 0.22217211134397474, 1.412231486567074, 1.6702402493281114, 0.1500938705500959],
        Precip        = [0.0, 0.0, 0.0, 0.9, 0.0, 7.7],
        Td            = [0.0, -2.0816681711721685e-17, -8.326672684688674e-17, 1.1848881316753612, 2.220446049250313e-16, 1.3877787807814457e-16],
        D             = [-0.22315088734244626, -0.21780089386670384, -2.1067555203802275, -0.07139512579375625, -1.0321020339500613, -0.2620446279917528],
        R             = [-0.0, -0.0, -0.0, -0.0, -0.0, -0.0],
        Swat          = [83.86682880676413, 83.8768752816219, 104.06567612548068, 55.97381850215723, 91.47999383954011, 76.46326291343296])
    # @test myDFcompare(df_partitioning_daily[[10, 20, 100, 180, 250, 260], cols_to_check], ref_partitioning_daily)
    @test isapprox(ref_partitioning_daily, df_partitioning_daily[[10, 20, 100, 180, 250, 260], cols_to_check], atol = 1e-5, rtol = 1e-5)

    # df_partitioning_monthly[:,cols_to_check]
    # julian_dput(df_partitioning_monthly[:,cols_to_check], multiline=true)
    ref_partitioning_monthly = DataFrame(
        ETa           = [12.608407286459945, 13.030233464936641, 18.324880172911683, 22.793362511281824, 44.55373047038666, 67.09866752231822, 65.25632926619858, 58.3262345696966, 47.177367338919545, 22.429812241895434, 10.87509771724212, 9.94579920104753, 0.1706586459056163],
        Esoil         = [0.19033084822860155, 0.0, 0.0, 1.7979749922902446, 5.518360669247835, 2.1357910332788848, 6.207931710760416, 6.581155515804109, 4.38340521073354, 3.0191932456006834, 0.8234998620788185, 0.0, 0.0],
        Esnow         = [1.0230481932789828, 3.5671213461847207, 5.699501531769204, 3.9862166438353714, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8944827966946712, 2.2122171998658007, 0.04886570426857504],
        Einterception = [10.669745233799594, 4.969997016210468, 6.735275083355589, 4.9930110467008575, 15.392808098042321, 10.556699999999998, 22.518735020696223, 17.09685796630591, 8.777594439949704, 2.3913881845359146, 5.507421255546405, 6.581318890654805, 0.0],
        Ta            = [0.7252830111527682, 4.493115102541451, 5.89010355778689, 12.016159828455352, 23.642561703096494, 54.40617648903933, 36.52966253474194, 34.64822108758658, 34.01636768823631, 17.01923081175883, 3.6496938029222235, 1.1522631105269225, 0.12179294163704127],
        Precip        = [158.1, 24.5, 50.4, 30.200000000000006, 91.30000000000001, 45.699999999999996, 159.40000000000003, 181.29999999999998, 56.20000000000001, 12.099999999999998, 83.7, 63.7, 0.0],
        Td            = [-7.003946034256359e-17, -3.5561831257524545e-16, -5.169475958410885e-16, -7.632783294297951e-16, 5.412337245047638e-16, 1.9743667924299457, 1.7208456881689926e-15, 1.457167719820518e-15, 7.494005416219807e-16, -1.4849232954361469e-15, 1.7520707107365752e-16, -2.1380466841414147e-16, 1.3877787807814457e-17],
        D             = [-6.5373064606680815, -44.23791134650466, -50.74748889094414, -90.21580154788323, -45.920744788272074, -8.952106091154258, -41.8123729912082, -129.34027158975402, -18.237643562579816, -14.12413925158462, -50.92796521872846, -12.760853999289004, -0.8137251051499813],
        R             = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        Swat          = [84.82260498498727, 102.99984953850444, 101.60690267140986, 104.71208507154076, 101.61283070822769, 72.0243763474135, 92.13562023506732, 106.7414216684709, 88.40709096931324, 86.57849430419505, 100.42812962836246, 89.50517207969122, 113.95770256535955])
    # @test myDFcompare(df_partitioning_monthly[:,cols_to_check], ref_partitioning_monthly)
    @test isapprox(ref_partitioning_monthly, df_partitioning_monthly[:,cols_to_check], atol = 1e-5, rtol = 1e-5)
    
    # julian_dput(df_partitioning_yearly[:,cols_to_check], multiline=true)
    ref_partitioning_yearly = DataFrame(
        ETa           = [392.41992176329455, 0.1706586459056163],
        Esoil         = [30.657643088023104, 0.0],
        Esnow         = [17.382587711628748, 0.04886570426857504],
        Einterception = [116.19085223579783, 0.0],
        Ta            = [228.18883872784517, 0.12179294163704127],
        Precip        = [956.5999999999999, 0.0],
        Td            = [1.9743667924299473, 1.3877787807814457e-17],
        D             = [-513.8146057385712, -0.8137251051499813],
        R             = [0.0, 0.0],
        Swat          = [94.258193681618, 113.95770256535955])
    # @test myDFcompare(df_partitioning_yearly[:,cols_to_check], ref_partitioning_yearly)
    @test isapprox(ref_partitioning_yearly, df_partitioning_yearly[:,cols_to_check], atol = 1e-5, rtol = 1e-5)


    # B) check consistency of undocumented (internal) output functions: 
    # Check soil depth indices:
    depths_to_test_mm = [100, 1000, 1200, 200, 300, 400, 150, ] # test unsorted input
    @test_logs (:warn, r"below simulation domain") LWFBrook90.intern___get_soil_idx(simulation, depths_to_test_mm)
    idx_to_read_out = LWFBrook90.intern___get_soil_idx(simulation, depths_to_test_mm)
    @test idx_to_read_out == Dict(
        150  => 3, 100  => 2, 200  => 4, 300  => 6, 400  => 8, 1000 => 20,
        1200 => 0) # 1200 is below the simulation domain!
    valid_idx_to_read_out = LWFBrook90.intern___get_soil_idx(simulation, depths_to_test_mm; only_valid_idxs = true)

    # Check isotope output for plotting:
    df_isotopePlot, RWUPlotlabel = LWFBrook90.intern___get_data_for_isotopePlot(simulation)
    @test RWUPlotlabel == "mean RWU depth\n(based on uptake only)"
    # julian_dput(df_isotopePlot[[1, 10, 100, 300], :], multiline=true)
    reference_isotopePlotPermutedDims_new = DataFrame(
        days_to_read_out_d    = [0.0, 9.0, 99.0, 299.0],
        col_PREC_amt_dense    = [0.0, 0.0, 0.0, 0.0],
        col_PREC_d18O_dense   = [-15.04, -15.04, -14.42, -22.6],
        col_PREC_d2H_dense    = [-111.96, -111.96, -102.2, -171.53],
        col_PREC_d18O         = [-15.04, -15.04, -14.42, -22.6],
        col_PREC_d2H          = [-111.96, -111.96, -102.2, -171.53],
        col_IRRIG_amt_dense   = [0.0, 0.0, 0.0, 0.0],
        col_IRRIG_d18O_dense  = [nothing, nothing, nothing, nothing],
        col_IRRIG_d2H_dense   = [nothing, nothing, nothing, nothing],
        col_IRRIG_d18O        = [nothing, nothing, nothing, nothing],
        col_IRRIG_d2H         = [nothing, nothing, nothing, nothing],
        col_INTS_d18O         = [-13.0, -15.04, -14.42, -22.6],
        col_INTR_d18O         = [-13.0, -15.04, -14.42, -22.6],
        col_SNOW_d18O         = [-13.0, -15.04, -14.42, -22.6],
        col_GWAT_d18O         = [-13.0, -10.609721281737363, -14.068076417219263, -9.128574715172391],
        col_RWU_d18O          = [-10.11111, -10.267518779484064, -14.349325713635785, -10.44957431491067],
        col_XYL_d18O          = [-10.11111, -10.111729973656493, -12.285771886666092, -10.23433689302966],
        col_INTS_d2H          = [-95.0, -111.96, -102.2, -171.53],
        col_INTR_d2H          = [-95.0, -111.96, -102.2, -171.53],
        col_SNOW_d2H          = [-95.0, -111.96, -102.2, -171.53],
        col_GWAT_d2H          = [-95.0, -91.78230625739805, -100.72512050987932, -60.41029202484412],
        col_RWU_d2H           = [-91.1111, -91.77261891033439, -101.63731872603199, -72.38068669347486],
        col_XYL_d2H           = [-91.1111, -91.11372225279625, -97.4623386184832, -70.0574265377251],
        col_RWU_centroid_mm   = [NaN, 43.22335331429762, 43.219875896063, 50.39351555017721],
        d18O_25               = [-10.11111, -10.299309982029751, -14.358553491579178, -10.52479258595642],
        d18O_75               = [-10.11111, -10.18764620901989, -14.333820176017783, -10.350912837250483],
        d18O_125              = [-10.11111, -10.164815722760252, -14.30783859661345, -10.247284958428379],
        d18O_175              = [-10.11111, -10.154855220780153, -14.28115307173558, -10.168678088798309],
        d18O_225              = [-10.11111, -10.148719732283729, -14.254324205726398, -10.097325814419253],
        d18O_275              = [-10.11111, -10.14440350207366, -14.228211999165085, -10.030066580188931],
        d18O_325              = [-10.11111, -10.14105490741634, -14.20416373453515, -9.965601616159162],
        d18O_375              = [-10.11111, -10.136781821464371, -14.184416970870567, -9.901303038141627],
        d18O_425              = [-10.11111, -10.120939691621501, -14.171744959643133, -9.830378647691756],
        d18O_475              = [-10.11111, -10.113523533724965, -14.159353984805149, -9.765037079436873],
        d18O_525              = [-10.11111, -10.111547533103538, -14.147101119695023, -9.698887794174317],
        d18O_575              = [-10.11111, -10.111173048474921, -14.13510255299859, -9.632291364565909],
        d18O_625              = [-10.11111, -10.11111753708149, -14.123572339431751, -9.566187059475617],
        d18O_675              = [-10.11111, -10.11111076779047, -14.112730513515107, -9.501355164842414],
        d18O_725              = [-10.11111, -10.111110067850742, -14.10277612940147, -9.438164994186687],
        d18O_775              = [-10.11111, -10.111110005266617, -14.093888647777565, -9.377505496118438],
        d18O_825              = [-10.11111, -10.111110000362245, -14.086230404924558, -9.320745986979693],
        d18O_875              = [-10.11111, -10.111110000022224, -14.079936082172184, -9.269182243739014],
        d18O_925              = [-10.11111, -10.111110000001208, -14.07508523675152, -9.223970492334875],
        d18O_975              = [-10.11111, -10.111110000000059, -14.07166293346204, -9.18624803768401],
        d18O_1025             = [-10.11111, -10.111110000000021, -14.069530303990993, -9.157556704179177],
        d18O_1075             = [-10.11111, -10.11111, -14.068470479481723, -9.140666111708875],
        d2H_25                = [-91.1111, -91.90703430174163, -101.70746656448462, -73.13280666627725],
        d2H_75                = [-91.1111, -91.43494480348973, -101.51694374449022, -71.40823744049233],
        d2H_125               = [-91.1111, -91.33828770056635, -101.32385052517911, -70.33002871193949],
        d2H_175               = [-91.1111, -91.29614325633239, -101.13355313474395, -69.53396080471373],
        d2H_225               = [-91.1111, -91.27018804525625, -100.95105036875468, -68.83431270410216],
        d2H_275               = [-91.1111, -91.25192969443353, -100.78234174809523, -68.19092494321217],
        d2H_325               = [-91.1111, -91.23776276947082, -100.63488122911426, -67.58303995105288],
        d2H_375               = [-91.1111, -91.21967134875656, -100.51961807890355, -66.98085892031288],
        d2H_425               = [-91.1111, -91.15271062199409, -100.44915859089882, -66.32180473627982],
        d2H_475               = [-91.1111, -91.12132932139045, -100.38588727419528, -65.71500685711791],
        d2H_525               = [-91.1111, -91.11295684379422, -100.32976311065804, -65.1034856181848],
        d2H_575               = [-91.1111, -91.11136793586604, -100.28270052788567, -64.49468963046915],
        d2H_625               = [-91.1111, -91.11113207491184, -100.24740741353222, -63.900786742744295],
        d2H_675               = [-91.1111, -91.11110327207437, -100.22668768445845, -63.33061959669466],
        d2H_725               = [-91.1111, -91.11110028957786, -100.223184249071, -62.78815078458296],
        d2H_775               = [-91.1111, -91.11110002251054, -100.23926954795631, -62.28129837008406],
        d2H_825               = [-91.1111, -91.11110000155067, -100.27682832706543, -61.82127780827766],
        d2H_875               = [-91.1111, -91.11110000009526, -100.33671144984096, -61.41753685382203],
        d2H_925               = [-91.1111, -91.11110000000522, -100.41759595196665, -61.076877076011854],
        d2H_975               = [-91.1111, -91.11110000000022, -100.51384971056115, -60.80393014634156],
        d2H_1025              = [-91.1111, -91.1111, -100.61184873934504, -60.60400050868431],
        d2H_1075              = [-91.1111, -91.1111, -100.68439413132077, -60.48945504507659]
    )
    reference_isotopePlotPermutedDims = permutedims(Matrix(reference_isotopePlotPermutedDims_new))    
    #@test myDFcompare(df_isotopePlot[[1, 10, 100, 300], :], reference_isotopePlotPermutedDims_new)
    @testset "postprocess_isotope-plot" begin
        for (it,(k,v)) in enumerate(pairs(eachcol(df_isotopePlot[[1, 10, 100, 300], :])))
            @test filter(!isnan, replace(v, nothing => NaN)) ≈ filter(!isnan, replace(reference_isotopePlotPermutedDims, nothing => NaN)[it,:]) ||
                "For $k: compared (now:) $v ≈ (reference:) $(replace(reference_isotopePlotPermutedDims, nothing => NaN)[it,:])"
        end
    end
end


@testset "simulate-and-postprocess" begin
    Δz_m = fill(0.05, 22)
    parametrizedSPAC = loadSPAC(
            "../examples/DAV2020-full/", "DAV2020-full"; 
            simulate_isotopes = true, simulate_irrigation = false,
            Δz_thickness_m = fill(0.05, 22),
            root_distribution = (beta = 0.77, z_rootMax_m = -0.5),
            IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -10.11111, delta2H_init_permil = -91.1111),
            canopy_evolution = (
                DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel = 100,
                LAI_rel = (DOY_Bstart = 120,
                           Bduration  = 21,
                           DOY_Cstart = 270,
                           Cduration  = 60,
                           LAI_perc_BtoC = 95,
                           LAI_perc_CtoB = 70)));

    simulation = setup(parametrizedSPAC);
    simulate!(simulation)

    # Check requirement of daily output resolution for some outputs:
    @test_throws r"Solution is not computed with daily output resolution." get_water_partitioning(simulation)
    # TODO: do we need that check for get_states, get_fluxes as well??

    # Rerun with daily:
    simulate!(simulation; 
        save_everystep = false, 
        saveat = range(parametrizedSPAC.tspan...), 
        tspan = parametrizedSPAC.tspan);

    # A) check consistency of documented output functions: 
    # (get_states, get_fluxes, get_forcing, get_soil_, get_water_partitioning)
    depths_to_test_mm_noWarning = [100, 1000, 200, 300, 400, 150, ]

    # compute current output to test
    df_get_forcing = get_forcing(simulation)[[1],:] # must be a DataFrame not a DataFrameRow
    df_get_states  = get_states(simulation, days_to_read_out_d = 0)
    df_get_fluxes  = get_fluxes(simulation, days_to_read_out_d = 0)
    df_get_soil_all= get_soil_([:θ,:ψ,:W,:SWATI,:K,:δ18O,:TRANI],
                               simulation, depths_to_read_out_mm = depths_to_test_mm_noWarning, days_to_read_out_d = 0)
    df_partitioning_daily, df_partitioning_monthly, df_partitioning_yearly = 
        get_water_partitioning(simulation)

    # get reference values of DataFrames (not DataFrameRows !) with: 
    function julian_dput(df::DataFrame; multiline = false)
        sep = ifelse(multiline, "\n","")
        cols = [repr(df[!, col]) for col in names(df)]
        code = "DataFrame($sep$(join(["$c=$v" for (c,v) in zip(names(df), cols)], ",$sep")))"
        println(code)
        return nothing
    end

    # julian_dput(df_get_forcing, multiline=true)
    ref_get_forcing  = DataFrame(
        dates=[DateTime("2021-01-01T00:00:00")],
        globrad_MJDayM2=[5.53],
        tmax_degC=[0.9],
        tmin_degC=[-10.1],
        vappres_kPa=[0.26],
        windspeed_ms=[1.5],
        prec_mmDay=[0.0],
        precdelta18O_permil=[-15.04],
        precdelta2H_permil=[-111.96],
        irrig_mmDay=[0.0],
        irrigdelta18O_permil=[nothing],
        irrigdelta2H_permil=[nothing],
        densef_percent=[100.0],
        height_m=[25.0],
        lai_m2m2=[2.0999999999999996],
        sai_m2m2=[1.0])

    # julian_dput(df_get_states, multiline=true)
    ref_get_states   = DataFrame(
        dates=[DateTime("2021-01-01T00:00:00")],
        GWAT_mm=[1.0],
        INTS_mm=[0.0],
        INTR_mm=[0.0],
        SNOW_mm=[0.0],
        CC_MJm2=[0.0],
        SNOWLQ_mm=[0.0],
        SWAT_mm=[85.24942394707108],
        GWAT_d18O=[-13.0],
        INTS_d18O=[-13.0],
        INTR_d18O=[-13.0],
        SNOW_d18O=[-13.0],
        XYLEM_d18O=[-10.11111],
        GWAT_d2H=[-95.0],
        INTS_d2H=[-95.0],
        INTR_d2H=[-95.0],
        SNOW_d2H=[-95.0],
        XYLEM_d2H=[-91.1111],
        d18O_permil_50mm=[-10.11111],
        d18O_permil_100mm=[-10.11111],
        d18O_permil_150mm=[-10.11111],
        d18O_permil_200mm=[-10.11111],
        d18O_permil_250mm=[-10.11111],
        d18O_permil_300mm=[-10.11111],
        d18O_permil_350mm=[-10.11111],
        d18O_permil_400mm=[-10.11111],
        d18O_permil_450mm=[-10.11111],
        d18O_permil_500mm=[-10.11111],
        d18O_permil_550mm=[-10.11111],
        d18O_permil_600mm=[-10.11111],
        d18O_permil_650mm=[-10.11111],
        d18O_permil_700mm=[-10.11111],
        d18O_permil_750mm=[-10.11111],
        d18O_permil_800mm=[-10.11111],
        d18O_permil_850mm=[-10.11111],
        d18O_permil_900mm=[-10.11111],
        d18O_permil_950mm=[-10.11111],
        d18O_permil_1000mm=[-10.11111],
        d18O_permil_1050mm=[-10.11111],
        d18O_permil_1100mm=[-10.11111],
        d2H_permil_50mm=[-91.1111],
        d2H_permil_100mm=[-91.1111],
        d2H_permil_150mm=[-91.1111],
        d2H_permil_200mm=[-91.1111],
        d2H_permil_250mm=[-91.1111],
        d2H_permil_300mm=[-91.1111],
        d2H_permil_350mm=[-91.1111],
        d2H_permil_400mm=[-91.1111],
        d2H_permil_450mm=[-91.1111],
        d2H_permil_500mm=[-91.1111],
        d2H_permil_550mm=[-91.1111],
        d2H_permil_600mm=[-91.1111],
        d2H_permil_650mm=[-91.1111],
        d2H_permil_700mm=[-91.1111],
        d2H_permil_750mm=[-91.1111],
        d2H_permil_800mm=[-91.1111],
        d2H_permil_850mm=[-91.1111],
        d2H_permil_900mm=[-91.1111],
        d2H_permil_950mm=[-91.1111],
        d2H_permil_1000mm=[-91.1111],
        d2H_permil_1050mm=[-91.1111],
        d2H_permil_1100mm=[-91.1111],
        θ_m3m3_50mm=[0.20058687987546137],
        θ_m3m3_100mm=[0.20058687987546137],
        θ_m3m3_150mm=[0.20058687987546137],
        θ_m3m3_200mm=[0.20058687987546137],
        θ_m3m3_250mm=[0.20058687987546137],
        θ_m3m3_300mm=[0.20058687987546137],
        θ_m3m3_350mm=[0.20058687987546137],
        θ_m3m3_400mm=[0.20058687987546137],
        θ_m3m3_450mm=[0.20058687987546137],
        θ_m3m3_500mm=[0.20058687987546137],
        θ_m3m3_550mm=[0.20058687987546137],
        θ_m3m3_600mm=[0.20058687987546137],
        θ_m3m3_650mm=[0.20058687987546137],
        θ_m3m3_700mm=[0.20058687987546137],
        θ_m3m3_750mm=[0.20058687987546137],
        θ_m3m3_800mm=[0.20058687987546137],
        θ_m3m3_850mm=[0.20058687987546137],
        θ_m3m3_900mm=[0.20058687987546137],
        θ_m3m3_950mm=[0.20058687987546137],
        θ_m3m3_1000mm=[0.2005868798754614],
        θ_m3m3_1050mm=[0.20058687987546137],
        θ_m3m3_1100mm=[0.20058687987546137],
        ψ_kPa_50mm=[-7.0000000000000036],
        ψ_kPa_100mm=[-7.0000000000000036],
        ψ_kPa_150mm=[-7.0000000000000036],
        ψ_kPa_200mm=[-7.0000000000000036],
        ψ_kPa_250mm=[-7.0000000000000036],
        ψ_kPa_300mm=[-7.0000000000000036],
        ψ_kPa_350mm=[-7.0000000000000036],
        ψ_kPa_400mm=[-7.0000000000000036],
        ψ_kPa_450mm=[-7.0000000000000036],
        ψ_kPa_500mm=[-7.0000000000000036],
        ψ_kPa_550mm=[-7.0000000000000036],
        ψ_kPa_600mm=[-7.0000000000000036],
        ψ_kPa_650mm=[-7.0000000000000036],
        ψ_kPa_700mm=[-7.0000000000000036],
        ψ_kPa_750mm=[-7.0000000000000036],
        ψ_kPa_800mm=[-7.0000000000000036],
        ψ_kPa_850mm=[-7.0000000000000036],
        ψ_kPa_900mm=[-7.0000000000000036],
        ψ_kPa_950mm=[-7.0000000000000036],
        ψ_kPa_1000mm=[-6.9999999999999964],
        ψ_kPa_1050mm=[-7.0000000000000036],
        ψ_kPa_1100mm=[-7.0000000000000036])
    
    # julian_dput(df_get_fluxes, multiline=true)
    ref_get_fluxes   = DataFrame(
        dates=[DateTime("2021-01-01T00:00:00")],
        cum_d_prec=[0.0],
        cum_d_sfal=[0.0],
        cum_d_sthr=[0.0],
        cum_d_sint=[0.0],
        cum_d_irrig=[0.0],
        cum_d_rfal=[0.0],
        cum_d_rint=[0.0],
        cum_d_rthr=[0.0],
        cum_d_rsno=[0.0],
        cum_d_rnet=[0.0],
        cum_d_smlt=[0.0],
        cum_d_irvp=[0.0],
        cum_d_isvp=[0.0],
        cum_d_snvp=[0.0],
        cum_d_slvp=[0.04482920026138777],
        cum_d_tran=[0.056616161386777504],
        RWU=[0.056616161386777504],
        INTS=[0.0],
        INTR=[0.0],
        SNOW=[0.0],
        CC=[0.0],
        SNOWLQ=[0.0],
        RWU_d18O=[-10.11111],
        RWU_d2H=[-91.1111],
        PREC_d18O=[-15.04],
        PREC_d2H=[-111.96],
        cum_d_pint=[2.6465243861228247],
        cum_d_ptran=[0.056616161386777504],
        cum_d_pslvp=[0.36472665495887147],
        srfl=[0.0],
        slfl=[0.0],
        byfl=[0.0],
        dsfl=[0.0],
        gwfl=[0.15412189586625977],
        vrfln=[0.15412189586625977],
        flow=[0.15412189586625977],
        seep=[0.0],
        evap=[0.10144536164816528],
        TRANI_mmday_50mm=[0.043498767636911864],
        TRANI_mmday_100mm=[0.009959724394658772],
        TRANI_mmday_150mm=[0.002409761574167726],
        TRANI_mmday_200mm=[0.0005748250905928008],
        TRANI_mmday_250mm=[0.0001346303525887302],
        TRANI_mmday_300mm=[3.0767571641898455e-5],
        TRANI_mmday_350mm=[6.792328768255243e-6],
        TRANI_mmday_400mm=[1.4228398921497332e-6],
        TRANI_mmday_450mm=[-5.304024446942134e-7],
        TRANI_mmday_500mm=[0.0],
        TRANI_mmday_550mm=[0.0],
        TRANI_mmday_600mm=[0.0],
        TRANI_mmday_650mm=[0.0],
        TRANI_mmday_700mm=[0.0],
        TRANI_mmday_750mm=[0.0],
        TRANI_mmday_800mm=[0.0],
        TRANI_mmday_850mm=[0.0],
        TRANI_mmday_900mm=[0.0],
        TRANI_mmday_950mm=[0.0],
        TRANI_mmday_1000mm=[0.0],
        TRANI_mmday_1050mm=[0.0],
        TRANI_mmday_1100mm=[0.0])
    
    # julian_dput(df_get_soil_all, multiline=true)
    ref_get_soil_all = DataFrame(
        time=[0], # TODO: note this is time as Int
        K_mmday_100mm=[1.1082021409141438],
        K_mmday_150mm=[1.1082021409141438],
        K_mmday_200mm=[1.1082021409141438],
        K_mmday_300mm=[1.1082021409141438],
        K_mmday_400mm=[1.1082021409141438],
        K_mmday_1000mm=[1.1082021409141432],
        SWATI_mm_100mm=[6.268339996108168],
        SWATI_mm_150mm=[6.268339996108169],
        SWATI_mm_200mm=[6.268339996108166],
        SWATI_mm_300mm=[6.268339996108166],
        SWATI_mm_400mm=[6.268339996108166],
        SWATI_mm_1000mm=[2.507335998443264],
        TRANI_mmday_100mm=[0.0],
        TRANI_mmday_150mm=[0.0],
        TRANI_mmday_200mm=[0.0],
        TRANI_mmday_300mm=[0.0],
        TRANI_mmday_400mm=[0.0],
        TRANI_mmday_1000mm=[0.0],
        W___100mm=[0.5298121496974679],
        W___150mm=[0.5298121496974679],
        W___200mm=[0.5298121496974679],
        W___300mm=[0.5298121496974679],
        W___400mm=[0.5298121496974679],
        W___1000mm=[0.529812149697468],
        δ18O_permil_100mm=[-10.11111],
        δ18O_permil_150mm=[-10.11111],
        δ18O_permil_200mm=[-10.11111],
        δ18O_permil_300mm=[-10.11111],
        δ18O_permil_400mm=[-10.11111],
        δ18O_permil_1000mm=[-10.11111],
        θ_m3m3_100mm=[0.20058687987546137],
        θ_m3m3_150mm=[0.20058687987546137],
        θ_m3m3_200mm=[0.20058687987546137],
        θ_m3m3_300mm=[0.20058687987546137],
        θ_m3m3_400mm=[0.20058687987546137],
        θ_m3m3_1000mm=[0.2005868798754614],
        ψ_kPa_100mm=[-7.0000000000000036],
        ψ_kPa_150mm=[-7.0000000000000036],
        ψ_kPa_200mm=[-7.0000000000000036],
        ψ_kPa_300mm=[-7.0000000000000036],
        ψ_kPa_400mm=[-7.0000000000000036],
        ψ_kPa_1000mm=[-6.9999999999999964])

    # compare current and reference output
    # set up approximate comparison of DataFrames (including nothgin, DateTime, ...)
    # https://discourse.julialang.org/t/is-there-a-package-to-compare-if-two-dataframes-are-the-same/44381/4
    mycompare(x::Nothing, y::Nothing; kwargs...) = true
    mycompare(x::Number, y::Number; rtol::Real=√eps(), atol::Real=0) = isapprox(x, y)
    mycompare(x, y; kwargs...) = isequal(x, y)
    function myDFcompare(x::DataFrame, y::DataFrame;
                      rtol::Real=√eps(), atol::Real=0) 
        print(x);print("\n")
        print(y);print("\n")
        all([all(vec) for vec in eachcol(mycompare.(x, y))])
    end

    @test myDFcompare(df_get_forcing, ref_get_forcing)
    @test myDFcompare(df_get_forcing, ref_get_forcing)
    @test myDFcompare(df_get_states,  ref_get_states)
    @test myDFcompare(df_get_fluxes,  ref_get_fluxes)
    # @test myDFcompare(df_get_soil_1,  ref_get_soil_1)
    # @test myDFcompare(df_get_soil_2,  ref_get_soil_2)
    # @test myDFcompare(df_get_soil_3,  ref_get_soil_3)
    # @test myDFcompare(df_get_soil_4,  ref_get_soil_4)
    # @test myDFcompare(df_get_soil_5,  ref_get_soil_5)
    # @test myDFcompare(df_get_soil_6,  ref_get_soil_6)
    # @test myDFcompare(df_get_soil_7,  ref_get_soil_7)

    # Check output of water partitioning output
    cols_to_check = ["ETa", "Esoil", "Esnow", "Einterception", "Ta", "Precip", "Td", "D", "R", "Swat",]
    # julian_dput(df_partitioning_daily[[10, 20, 100, 180, 250, 260], cols_to_check], multiline=true)
    ref_partitioning_daily = DataFrame(
        ETa           = [0.041923240612276605, 0.04854475293069882, 0.4374160550376799, 1.6240215401844946, 1.8509559669497675, 1.2653333329108203],
        Esoil         = [0.009285599165961049, 0.0, 0.0, 0.0038900536174205622, 0.18071571762165603, 0.027720901826861828],
        Esnow         = [0.0, 0.016611946626641515, 0.21524394369370514, 0.0, 0.0, 0.0],
        Einterception = [0.0, 0.0, 0.0, 0.2079, 0.0, 1.0875185605338626],
        Ta            = [0.032637641446315556, 0.0319328063040573, 0.22217211134397474, 1.412231486567074, 1.6702402493281114, 0.1500938705500959],
        Precip        = [0.0, 0.0, 0.0, 0.9, 0.0, 7.7],
        Td            = [0.0, -2.0816681711721685e-17, -8.326672684688674e-17, 1.1848881316753612, 2.220446049250313e-16, 1.3877787807814457e-16],
        D             = [-0.22315088734244626, -0.21780089386670384, -2.1067555203802275, -0.07139512579375625, -1.0321020339500613, -0.2620446279917528],
        R             = [-0.0, -0.0, -0.0, -0.0, -0.0, -0.0],
        Swat          = [83.86682880676413, 83.8768752816219, 104.06567612548068, 55.97381850215723, 91.47999383954011, 76.46326291343296])
    # @test myDFcompare(df_partitioning_daily[[10, 20, 100, 180, 250, 260], cols_to_check], ref_partitioning_daily)
    @test isapprox(ref_partitioning_daily, df_partitioning_daily[[10, 20, 100, 180, 250, 260], cols_to_check], atol = 1e-5, rtol = 1e-5)

    # df_partitioning_monthly[:,cols_to_check]
    # julian_dput(df_partitioning_monthly[:,cols_to_check], multiline=true)
    ref_partitioning_monthly = DataFrame(
        ETa           = [12.608407286459945, 13.030233464936641, 18.324880172911683, 22.793362511281824, 44.55373047038666, 67.09866752231822, 65.25632926619858, 58.3262345696966, 47.177367338919545, 22.429812241895434, 10.87509771724212, 9.94579920104753, 0.1706586459056163],
        Esoil         = [0.19033084822860155, 0.0, 0.0, 1.7979749922902446, 5.518360669247835, 2.1357910332788848, 6.207931710760416, 6.581155515804109, 4.38340521073354, 3.0191932456006834, 0.8234998620788185, 0.0, 0.0],
        Esnow         = [1.0230481932789828, 3.5671213461847207, 5.699501531769204, 3.9862166438353714, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8944827966946712, 2.2122171998658007, 0.04886570426857504],
        Einterception = [10.669745233799594, 4.969997016210468, 6.735275083355589, 4.9930110467008575, 15.392808098042321, 10.556699999999998, 22.518735020696223, 17.09685796630591, 8.777594439949704, 2.3913881845359146, 5.507421255546405, 6.581318890654805, 0.0],
        Ta            = [0.7252830111527682, 4.493115102541451, 5.89010355778689, 12.016159828455352, 23.642561703096494, 54.40617648903933, 36.52966253474194, 34.64822108758658, 34.01636768823631, 17.01923081175883, 3.6496938029222235, 1.1522631105269225, 0.12179294163704127],
        Precip        = [158.1, 24.5, 50.4, 30.200000000000006, 91.30000000000001, 45.699999999999996, 159.40000000000003, 181.29999999999998, 56.20000000000001, 12.099999999999998, 83.7, 63.7, 0.0],
        Td            = [-7.003946034256359e-17, -3.5561831257524545e-16, -5.169475958410885e-16, -7.632783294297951e-16, 5.412337245047638e-16, 1.9743667924299457, 1.7208456881689926e-15, 1.457167719820518e-15, 7.494005416219807e-16, -1.4849232954361469e-15, 1.7520707107365752e-16, -2.1380466841414147e-16, 1.3877787807814457e-17],
        D             = [-6.5373064606680815, -44.23791134650466, -50.74748889094414, -90.21580154788323, -45.920744788272074, -8.952106091154258, -41.8123729912082, -129.34027158975402, -18.237643562579816, -14.12413925158462, -50.92796521872846, -12.760853999289004, -0.8137251051499813],
        R             = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        Swat          = [84.82260498498727, 102.99984953850444, 101.60690267140986, 104.71208507154076, 101.61283070822769, 72.0243763474135, 92.13562023506732, 106.7414216684709, 88.40709096931324, 86.57849430419505, 100.42812962836246, 89.50517207969122, 113.95770256535955])
    # @test myDFcompare(df_partitioning_monthly[:,cols_to_check], ref_partitioning_monthly)
    @test isapprox(ref_partitioning_monthly, df_partitioning_monthly[:,cols_to_check], atol = 1e-5, rtol = 1e-5)
    
    # julian_dput(df_partitioning_yearly[:,cols_to_check], multiline=true)
    ref_partitioning_yearly = DataFrame(
        ETa           = [392.41992176329455, 0.1706586459056163],
        Esoil         = [30.657643088023104, 0.0],
        Esnow         = [17.382587711628748, 0.04886570426857504],
        Einterception = [116.19085223579783, 0.0],
        Ta            = [228.18883872784517, 0.12179294163704127],
        Precip        = [956.5999999999999, 0.0],
        Td            = [1.9743667924299473, 1.3877787807814457e-17],
        D             = [-513.8146057385712, -0.8137251051499813],
        R             = [0.0, 0.0],
        Swat          = [94.258193681618, 113.95770256535955])
    # @test myDFcompare(df_partitioning_yearly[:,cols_to_check], ref_partitioning_yearly)
    @test isapprox(ref_partitioning_yearly, df_partitioning_yearly[:,cols_to_check], atol = 1e-5, rtol = 1e-5)


    # B) check consistency of undocumented (internal) output functions: 
    # Check soil depth indices:
    depths_to_test_mm = [100, 1000, 1200, 200, 300, 400, 150, ] # test unsorted input
    @test_logs (:warn, r"below simulation domain") LWFBrook90.intern___get_soil_idx(simulation, depths_to_test_mm)
    idx_to_read_out = LWFBrook90.intern___get_soil_idx(simulation, depths_to_test_mm)
    @test idx_to_read_out == Dict(
        150  => 3, 100  => 2, 200  => 4, 300  => 6, 400  => 8, 1000 => 20,
        1200 => 0) # 1200 is below the simulation domain!
    valid_idx_to_read_out = LWFBrook90.intern___get_soil_idx(simulation, depths_to_test_mm; only_valid_idxs = true)

    # Check isotope output for plotting:
    df_isotopePlot, RWUPlotlabel = LWFBrook90.intern___get_data_for_isotopePlot(simulation)
    @test RWUPlotlabel == "mean RWU depth\n(based on uptake only)"
    # julian_dput(df_isotopePlot[[1, 10, 100, 300], :], multiline=true)
    reference_isotopePlotPermutedDims_new = DataFrame(
        days_to_read_out_d    = [0.0, 9.0, 99.0, 299.0],
        col_PREC_amt_dense    = [0.0, 0.0, 0.0, 0.0],
        col_PREC_d18O_dense   = [-15.04, -15.04, -14.42, -22.6],
        col_PREC_d2H_dense    = [-111.96, -111.96, -102.2, -171.53],
        col_PREC_d18O         = [-15.04, -15.04, -14.42, -22.6],
        col_PREC_d2H          = [-111.96, -111.96, -102.2, -171.53],
        col_IRRIG_amt_dense   = [0.0, 0.0, 0.0, 0.0],
        col_IRRIG_d18O_dense  = [nothing, nothing, nothing, nothing],
        col_IRRIG_d2H_dense   = [nothing, nothing, nothing, nothing],
        col_IRRIG_d18O        = [nothing, nothing, nothing, nothing],
        col_IRRIG_d2H         = [nothing, nothing, nothing, nothing],
        col_INTS_d18O         = [-13.0, -15.04, -14.42, -22.6],
        col_INTR_d18O         = [-13.0, -15.04, -14.42, -22.6],
        col_SNOW_d18O         = [-13.0, -15.04, -14.42, -22.6],
        col_GWAT_d18O         = [-13.0, -10.609721281737363, -14.068076417219263, -9.128574715172391],
        col_RWU_d18O          = [-10.11111, -10.267518779484064, -14.349325713635785, -10.44957431491067],
        col_XYL_d18O          = [-10.11111, -10.111729973656493, -12.285771886666092, -10.23433689302966],
        col_INTS_d2H          = [-95.0, -111.96, -102.2, -171.53],
        col_INTR_d2H          = [-95.0, -111.96, -102.2, -171.53],
        col_SNOW_d2H          = [-95.0, -111.96, -102.2, -171.53],
        col_GWAT_d2H          = [-95.0, -91.78230625739805, -100.72512050987932, -60.41029202484412],
        col_RWU_d2H           = [-91.1111, -91.77261891033439, -101.63731872603199, -72.38068669347486],
        col_XYL_d2H           = [-91.1111, -91.11372225279625, -97.4623386184832, -70.0574265377251],
        col_RWU_centroid_mm   = [NaN, 43.22335331429762, 43.219875896063, 50.39351555017721],
        d18O_25               = [-10.11111, -10.299309982029751, -14.358553491579178, -10.52479258595642],
        d18O_75               = [-10.11111, -10.18764620901989, -14.333820176017783, -10.350912837250483],
        d18O_125              = [-10.11111, -10.164815722760252, -14.30783859661345, -10.247284958428379],
        d18O_175              = [-10.11111, -10.154855220780153, -14.28115307173558, -10.168678088798309],
        d18O_225              = [-10.11111, -10.148719732283729, -14.254324205726398, -10.097325814419253],
        d18O_275              = [-10.11111, -10.14440350207366, -14.228211999165085, -10.030066580188931],
        d18O_325              = [-10.11111, -10.14105490741634, -14.20416373453515, -9.965601616159162],
        d18O_375              = [-10.11111, -10.136781821464371, -14.184416970870567, -9.901303038141627],
        d18O_425              = [-10.11111, -10.120939691621501, -14.171744959643133, -9.830378647691756],
        d18O_475              = [-10.11111, -10.113523533724965, -14.159353984805149, -9.765037079436873],
        d18O_525              = [-10.11111, -10.111547533103538, -14.147101119695023, -9.698887794174317],
        d18O_575              = [-10.11111, -10.111173048474921, -14.13510255299859, -9.632291364565909],
        d18O_625              = [-10.11111, -10.11111753708149, -14.123572339431751, -9.566187059475617],
        d18O_675              = [-10.11111, -10.11111076779047, -14.112730513515107, -9.501355164842414],
        d18O_725              = [-10.11111, -10.111110067850742, -14.10277612940147, -9.438164994186687],
        d18O_775              = [-10.11111, -10.111110005266617, -14.093888647777565, -9.377505496118438],
        d18O_825              = [-10.11111, -10.111110000362245, -14.086230404924558, -9.320745986979693],
        d18O_875              = [-10.11111, -10.111110000022224, -14.079936082172184, -9.269182243739014],
        d18O_925              = [-10.11111, -10.111110000001208, -14.07508523675152, -9.223970492334875],
        d18O_975              = [-10.11111, -10.111110000000059, -14.07166293346204, -9.18624803768401],
        d18O_1025             = [-10.11111, -10.111110000000021, -14.069530303990993, -9.157556704179177],
        d18O_1075             = [-10.11111, -10.11111, -14.068470479481723, -9.140666111708875],
        d2H_25                = [-91.1111, -91.90703430174163, -101.70746656448462, -73.13280666627725],
        d2H_75                = [-91.1111, -91.43494480348973, -101.51694374449022, -71.40823744049233],
        d2H_125               = [-91.1111, -91.33828770056635, -101.32385052517911, -70.33002871193949],
        d2H_175               = [-91.1111, -91.29614325633239, -101.13355313474395, -69.53396080471373],
        d2H_225               = [-91.1111, -91.27018804525625, -100.95105036875468, -68.83431270410216],
        d2H_275               = [-91.1111, -91.25192969443353, -100.78234174809523, -68.19092494321217],
        d2H_325               = [-91.1111, -91.23776276947082, -100.63488122911426, -67.58303995105288],
        d2H_375               = [-91.1111, -91.21967134875656, -100.51961807890355, -66.98085892031288],
        d2H_425               = [-91.1111, -91.15271062199409, -100.44915859089882, -66.32180473627982],
        d2H_475               = [-91.1111, -91.12132932139045, -100.38588727419528, -65.71500685711791],
        d2H_525               = [-91.1111, -91.11295684379422, -100.32976311065804, -65.1034856181848],
        d2H_575               = [-91.1111, -91.11136793586604, -100.28270052788567, -64.49468963046915],
        d2H_625               = [-91.1111, -91.11113207491184, -100.24740741353222, -63.900786742744295],
        d2H_675               = [-91.1111, -91.11110327207437, -100.22668768445845, -63.33061959669466],
        d2H_725               = [-91.1111, -91.11110028957786, -100.223184249071, -62.78815078458296],
        d2H_775               = [-91.1111, -91.11110002251054, -100.23926954795631, -62.28129837008406],
        d2H_825               = [-91.1111, -91.11110000155067, -100.27682832706543, -61.82127780827766],
        d2H_875               = [-91.1111, -91.11110000009526, -100.33671144984096, -61.41753685382203],
        d2H_925               = [-91.1111, -91.11110000000522, -100.41759595196665, -61.076877076011854],
        d2H_975               = [-91.1111, -91.11110000000022, -100.51384971056115, -60.80393014634156],
        d2H_1025              = [-91.1111, -91.1111, -100.61184873934504, -60.60400050868431],
        d2H_1075              = [-91.1111, -91.1111, -100.68439413132077, -60.48945504507659]
    )
    reference_isotopePlotPermutedDims = permutedims(Matrix(reference_isotopePlotPermutedDims_new))    
    #@test myDFcompare(df_isotopePlot[[1, 10, 100, 300], :], reference_isotopePlotPermutedDims_new)
    @testset "postprocess_isotope-plot" begin
        for (it,(k,v)) in enumerate(pairs(eachcol(df_isotopePlot[[1, 10, 100, 300], :])))
            @test filter(!isnan, replace(v, nothing => NaN)) ≈ filter(!isnan, replace(reference_isotopePlotPermutedDims, nothing => NaN)[it,:]) ||
                "For $k: compared (now:) $v ≈ (reference:) $(replace(reference_isotopePlotPermutedDims, nothing => NaN)[it,:])"
        end
    end
end
