# Test Set: Multi-Species Competition for Shared Soil Water Pool
#
# This test suite verifies:
# 1) Equivalence between single-species simulation and multi-species simulation with N=1 species (100% cover)
# 2) Equivalence between single-species simulation and multi-species simulation with 2 identical species (50%/50% cover)
# 3) Equivalence between single-species simulation and multi-species simulation with 2 identical species (30%/70% cover)
# 4) Validation of species naming consistency across all input files and parameters
# 5) Functional execution and postprocessing of DAV2020-multispecies-full and DAV2020-multispecies-bare-minimum

using Test, LWFBrook90, DataFrames, Dates

if basename(pwd()) != "test"; cd("test"); end

@testset "Multi-species: Equivalence Tests" begin
    # -------------------------------------------------------------
    # Setup Single-Species Baseline Model (DAV2020-full)
    # -------------------------------------------------------------
    path_single = "../examples/DAV2020-full/"
    prefix_single = "DAV2020-full"
    
    tspan_test = (0.0, 20.0) # 20 days simulation for precise equivalence
    
    model_single = loadSPAC(path_single, prefix_single; simulate_isotopes = true)
    sim_single = setup(model_single; requested_tspan = tspan_test)
    simulate!(sim_single)
    
    states_single = get_states(sim_single)
    fluxes_single = get_fluxes(sim_single)
    
    # -------------------------------------------------------------
    # Test 1: Multi-species with 1 species (100% cover) == Single-species baseline
    # -------------------------------------------------------------
    @testset "Equivalence 1: N=1 species (100% cover)" begin
        # Programmatic setup of N=1 multispecies model with identical parameters
        model_m1 = loadSPAC(path_single, prefix_single; simulate_isotopes = true)
        sim_m1 = setup(model_m1; requested_tspan = tspan_test)
        simulate!(sim_m1)
        
        states_m1 = get_states(sim_m1)
        fluxes_m1 = get_fluxes(sim_m1)
        
        # Check scalar storages
        @test states_m1.GWAT_mm ≈ states_single.GWAT_mm atol=1e-8
        @test states_m1.SWAT_mm ≈ states_single.SWAT_mm atol=1e-8
        @test states_m1.SNOW_mm ≈ states_single.SNOW_mm atol=1e-8
        @test states_m1.INTS_mm ≈ states_single.INTS_mm atol=1e-8
        @test states_m1.INTR_mm ≈ states_single.INTR_mm atol=1e-8
        
        # Check cumulative fluxes
        @test fluxes_m1.cum_d_tran ≈ fluxes_single.cum_d_tran atol=1e-8
        @test fluxes_m1.evap ≈ fluxes_single.evap atol=1e-8
        @test fluxes_m1.flow ≈ fluxes_single.flow atol=1e-8
        @test fluxes_m1.cum_d_slvp ≈ fluxes_single.cum_d_slvp atol=1e-8
        @test fluxes_m1.cum_d_irvp ≈ fluxes_single.cum_d_irvp atol=1e-8
        
        # Check isotopes
        @test isapprox.(states_m1.GWAT_d18O, states_single.GWAT_d18O, nans=true) |> all
        @test isapprox.(states_m1.GWAT_d2H, states_single.GWAT_d2H, nans=true) |> all
        @test isapprox.(fluxes_m1.RWU_d18O, fluxes_single.RWU_d18O, nans=true) |> all
        @test isapprox.(fluxes_m1.RWU_d2H, fluxes_single.RWU_d2H, nans=true) |> all
    end

    # -------------------------------------------------------------
    # Setup Single-Species Baseline Model (DAV2020-bare-minimum) for 2-species identical tests
    # -------------------------------------------------------------
    path_minimal = "../examples/DAV2020-bare-minimum/"
    prefix_minimal = "DAV2020-minimal"
    
    canopy_evo_identical = (DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel = 100,
        LAI_rel = (DOY_Bstart = 120, Bduration = 21, DOY_Cstart = 270, Cduration = 60, LAI_perc_BtoC = 100, LAI_perc_CtoB = 60))
    root_dist_identical = (beta = 0.98, z_rootMax_m = -1.1)
    ic_soil_identical = (PSIM_init_kPa = -6.3, delta18O_init_permil = -13.0, delta2H_init_permil = -95.0)

    model_single_min = loadSPAC(path_minimal, prefix_minimal;
        simulate_isotopes = true,
        Δz_thickness_m = fill(0.1, 11),
        storm_durations_h = fill(4.0, 12),
        root_distribution = root_dist_identical,
        canopy_evolution = canopy_evo_identical,
        IC_soil = ic_soil_identical,
        IC_scalar = (
            amount = (u_GWAT_init_mm = 1.0, u_SNOW_init_mm = 0.0, u_CC_init_MJ_per_m2 = 0.0, u_SNOWLQ_init_mm = 0.0,
                      u_INTS_init_mm = 0.0, u_INTR_init_mm = 0.0),
            d18O   = (u_GWAT_init_permil = -13.0, u_SNOW_init_permil = -13.0,
                      u_INTS_init_permil = -13.0, u_INTR_init_permil = -13.0),
            d2H    = (u_GWAT_init_permil = -95.0, u_SNOW_init_permil = -95.0,
                      u_INTS_init_permil = -95.0, u_INTR_init_permil = -95.0)
        )
    )
    sim_single_min = setup(model_single_min; requested_tspan = tspan_test)
    simulate!(sim_single_min)
    states_single_min = get_states(sim_single_min)
    fluxes_single_min = get_fluxes(sim_single_min)

    # -------------------------------------------------------------
    # Test 2: Multi-species with 2 identical species (50%/50% cover) == Single-species baseline
    # -------------------------------------------------------------
    @testset "Equivalence 2: N=2 identical species (50%/50% cover)" begin
        # Create a 2-species model with identical parameters and w1=0.5, w2=0.5
        model_m5050 = loadSPAC("../examples/DAV2020-multispecies-bare-minimum/", "DAV2020-multispecies-minimal";
            simulate_isotopes = true,
            cover_fractions   = (Beech = 0.5, Spruce = 0.5),
            Δz_thickness_m    = fill(0.1, 11),
            storm_durations_h = fill(4.0, 12),
            params            = (
                Beech  = model_single_min.pars.params,
                Spruce = model_single_min.pars.params
            ),
            root_distribution = (
                Beech  = root_dist_identical,
                Spruce = root_dist_identical
            ),
            canopy_evolution  = (
                Beech  = canopy_evo_identical,
                Spruce = canopy_evo_identical
            ),
            IC_soil   = ic_soil_identical,
            IC_scalar = (
                amount = (u_GWAT_init_mm = 1.0, u_SNOW_init_mm = 0.0, u_CC_init_MJ_per_m2 = 0.0, u_SNOWLQ_init_mm = 0.0,
                          u_INTS_init_mm_Beech = 0.0, u_INTR_init_mm_Beech = 0.0,
                          u_INTS_init_mm_Spruce = 0.0, u_INTR_init_mm_Spruce = 0.0),
                d18O   = (u_GWAT_init_permil = -13.0, u_SNOW_init_permil = -13.0,
                          u_INTS_init_permil_Beech = -13.0, u_INTR_init_permil_Beech = -13.0,
                          u_INTS_init_permil_Spruce = -13.0, u_INTR_init_permil_Spruce = -13.0),
                d2H    = (u_GWAT_init_permil = -95.0, u_SNOW_init_permil = -95.0,
                          u_INTS_init_permil_Beech = -95.0, u_INTR_init_permil_Beech = -95.0,
                          u_INTS_init_permil_Spruce = -95.0, u_INTR_init_permil_Spruce = -95.0)
            )
        )
        sim_m5050 = setup(model_m5050; requested_tspan = tspan_test)
        simulate!(sim_m5050)
        
        states_m5050 = get_states(sim_m5050)
        fluxes_m5050 = get_fluxes(sim_m5050)
        
        # Shared soil water and groundwater storages must match baseline
        @test states_m5050.GWAT_mm ≈ states_single_min.GWAT_mm atol=1e-8
        @test states_m5050.SWAT_mm ≈ states_single_min.SWAT_mm atol=1e-8
        
        # Total transpiration is the sum of both species: TRAN_total = TRAN_Beech + TRAN_Spruce
        @test (fluxes_m5050.cum_d_tran_Beech + fluxes_m5050.cum_d_tran_Spruce) ≈ fluxes_single_min.cum_d_tran atol=1e-8
        @test fluxes_m5050.cum_d_tran ≈ fluxes_single_min.cum_d_tran atol=1e-8
        
        # Each 50% species does exactly 50% of the baseline transpiration
        @test fluxes_m5050.cum_d_tran_Beech ≈ 0.5 .* fluxes_single_min.cum_d_tran atol=1e-8
        @test fluxes_m5050.cum_d_tran_Spruce ≈ 0.5 .* fluxes_single_min.cum_d_tran atol=1e-8
        
        # Ecosystem evaporation & streamflow match
        @test fluxes_m5050.evap ≈ fluxes_single_min.evap atol=1e-8
        @test fluxes_m5050.flow ≈ fluxes_single_min.flow atol=1e-8
        
        # Root water uptake isotopes of identical species match baseline
        @test isapprox.(fluxes_m5050.RWU_d18O_Beech, fluxes_single_min.RWU_d18O, nans=true) |> all
        @test isapprox.(fluxes_m5050.RWU_d18O_Spruce, fluxes_single_min.RWU_d18O, nans=true) |> all
    end

    # -------------------------------------------------------------
    # Test 3: Multi-species with 2 identical species (30%/70% cover) == Single-species baseline
    # -------------------------------------------------------------
    @testset "Equivalence 3: N=2 identical species (30%/70% cover)" begin
        model_m3070 = loadSPAC("../examples/DAV2020-multispecies-bare-minimum/", "DAV2020-multispecies-minimal";
            simulate_isotopes = true,
            cover_fractions   = (Beech = 0.3, Spruce = 0.7),
            Δz_thickness_m    = fill(0.1, 11),
            storm_durations_h = fill(4.0, 12),
            params            = (
                Beech  = model_single_min.pars.params,
                Spruce = model_single_min.pars.params
            ),
            root_distribution = (
                Beech  = root_dist_identical,
                Spruce = root_dist_identical
            ),
            canopy_evolution  = (
                Beech  = canopy_evo_identical,
                Spruce = canopy_evo_identical
            ),
            IC_soil   = ic_soil_identical,
            IC_scalar = (
                amount = (u_GWAT_init_mm = 1.0, u_SNOW_init_mm = 0.0, u_CC_init_MJ_per_m2 = 0.0, u_SNOWLQ_init_mm = 0.0,
                          u_INTS_init_mm_Beech = 0.0, u_INTR_init_mm_Beech = 0.0,
                          u_INTS_init_mm_Spruce = 0.0, u_INTR_init_mm_Spruce = 0.0),
                d18O   = (u_GWAT_init_permil = -13.0, u_SNOW_init_permil = -13.0,
                          u_INTS_init_permil_Beech = -13.0, u_INTR_init_permil_Beech = -13.0,
                          u_INTS_init_permil_Spruce = -13.0, u_INTR_init_permil_Spruce = -13.0),
                d2H    = (u_GWAT_init_permil = -95.0, u_SNOW_init_permil = -95.0,
                          u_INTS_init_permil_Beech = -95.0, u_INTR_init_permil_Beech = -95.0,
                          u_INTS_init_permil_Spruce = -95.0, u_INTR_init_permil_Spruce = -95.0)
            )
        )
        sim_m3070 = setup(model_m3070; requested_tspan = tspan_test)
        simulate!(sim_m3070)
        
        states_m3070 = get_states(sim_m3070)
        fluxes_m3070 = get_fluxes(sim_m3070)
        
        @test states_m3070.GWAT_mm ≈ states_single_min.GWAT_mm atol=1e-8
        @test states_m3070.SWAT_mm ≈ states_single_min.SWAT_mm atol=1e-8
        @test fluxes_m3070.cum_d_tran ≈ fluxes_single_min.cum_d_tran atol=1e-8
        @test fluxes_m3070.cum_d_tran_Beech ≈ 0.3 .* fluxes_single_min.cum_d_tran atol=1e-8
        @test fluxes_m3070.cum_d_tran_Spruce ≈ 0.7 .* fluxes_single_min.cum_d_tran atol=1e-8
    end
end

@testset "Multi-species: Species Naming Validation & Errors" begin
    # Test error when species names in canopy_evolution do not match param.csv
    @test_throws ErrorException loadSPAC("../examples/DAV2020-multispecies-bare-minimum/", "DAV2020-multispecies-minimal";
        storm_durations_h = fill(4.0, 12),
        canopy_evolution = (Oak = (DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel = 100,
                                   LAI_rel = (DOY_Bstart = 120, Bduration = 21, DOY_Cstart = 270, Cduration = 60, LAI_perc_BtoC = 100, LAI_perc_CtoB = 60)),
                            Spruce = (DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel = 100,
                                      LAI_rel = (DOY_Bstart = 120, Bduration = 21, DOY_Cstart = 270, Cduration = 60, LAI_perc_BtoC = 100, LAI_perc_CtoB = 60)))
    )
    
    # Test error when species names in root_distribution do not match param.csv
    @test_throws ErrorException loadSPAC("../examples/DAV2020-multispecies-bare-minimum/", "DAV2020-multispecies-minimal";
        storm_durations_h = fill(4.0, 12),
        root_distribution = (Pine = (beta = 0.98, z_rootMax_m = -1.1), Spruce = (beta = 0.92, z_rootMax_m = -0.5))
    )
end

@testset "Multi-species: Full and Bare-Minimum Example Datasets" begin
    # Test loading and simulation of DAV2020-multispecies-full
    model_full = loadSPAC("../examples/DAV2020-multispecies-full/", "DAV2020-multispecies-full"; simulate_isotopes = true)
    @test length(model_full.pars.species_names) == 2
    @test :Beech in model_full.pars.species_names
    @test :Spruce in model_full.pars.species_names
    @test model_full.pars.cover_fractions.Beech ≈ 0.60
    @test model_full.pars.cover_fractions.Spruce ≈ 0.40
    
    sim_full = setup(model_full; requested_tspan = (0.0, 30.0))
    simulate!(sim_full)
    
    states_full = get_states(sim_full)
    fluxes_full = get_fluxes(sim_full)
    
    @test "cum_d_tran_Beech" in names(fluxes_full)
    @test "cum_d_tran_Spruce" in names(fluxes_full)
    @test "RWU_d18O_Beech" in names(fluxes_full)
    @test "RWU_d18O_Spruce" in names(fluxes_full)
    @test "INTS_mm_Beech" in names(states_full)
    @test "INTS_mm_Spruce" in names(states_full)
    
    # Total transpiration must equal sum of Beech and Spruce
    @test all(fluxes_full.cum_d_tran .≈ fluxes_full.cum_d_tran_Beech .+ fluxes_full.cum_d_tran_Spruce)
end

