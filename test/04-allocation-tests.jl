# Small hot-path regression checks; also runnable independently of plotting tests.
function norm_allocated(u)
    LWFBrook90.norm_to_use(u, 0.0) # compile before measuring
    return @allocated LWFBrook90.norm_to_use(u, 0.0)
end

function isotope_callback_allocated(integrator)
    callback = LWFBrook90.LWFBrook90R_updateIsotopes_GWAT_SWAT_AdvecDiff!
    callback(integrator.u, integrator.t, integrator) # compile before measuring
    return @allocated callback(integrator.u, integrator.t, integrator)
end

@testset "Per-step isotope callback: reuse caches and preserve root outfluxes" begin
    for fractionates in (false, true)
        model = loadSPAC(joinpath(@__DIR__, "..", "examples", "DAV2020-full"),
            "DAV2020-full"; simulate_isotopes=true,
            simulate_evaporation_fractionation=fractionates)
        simulation = remakeSPAC(model; params=(SLVPDEPTH_m=0.1,),
            requested_tspan=(0.0, 3.0))
        simulate!(simulation; progress=false, save_everystep=false, saveat=0.0:1.0:3.0)
        @test SciMLBase.successful_retcode(simulation.ODESolution)

        for uptake in (:positive, :negative, :zero), groundwater_mm in (0.0, 1.0)
            u = copy(simulation.ODESolution.u[end])
            p = deepcopy(simulation.ODEProblem.p) # never share scratch caches between runs
            u.GWAT.mm = groundwater_mm
            u.GWAT.d18O = -12.0
            u.GWAT.d2H = -85.0
            u.SWATI.d18O .= range(-15.0, -5.0; length=p.NLAYER)
            u.SWATI.d2H .= range(-100.0, -50.0; length=p.NLAYER)
            p.aux_du_TRANI .= uptake == :zero ? 0.0 : 0.001
            uptake == :negative && (p.aux_du_TRANI[1] = -0.0005)
            uprev = copy(u)
            uprev.SWATI.mm .*= 0.999
            integrator = (; u, uprev, p, t=3.0, tprev=3.0-1/240)
            θ_expected = LWFBrook90.KPT.derive_auxiliary_SOILVAR(u.SWATI.mm, p.p_soil)[4]
            xylem_before = (u.XYLEM.d18O, u.XYLEM.d2H)
            expected_rwu = map(((δsoil, δxylem, ratio),) -> begin
                concentrations = ifelse.(p.aux_du_TRANI .< 0,
                    LWFBrook90.ISO.δ_to_x(δxylem, ratio),
                    LWFBrook90.ISO.δ_to_x.(δsoil, ratio))
                # Original implementation, including signed weights and zero uptake.
                LWFBrook90.ISO.x_to_δ(
                    LWFBrook90.mean(concentrations, LWFBrook90.weights(p.aux_du_TRANI)), ratio)
            end, ((uprev.SWATI.d18O, u.XYLEM.d18O, LWFBrook90.ISO.R_VSMOW¹⁸O),
                  (uprev.SWATI.d2H, u.XYLEM.d2H, LWFBrook90.ISO.R_VSMOW²H)))

            LWFBrook90.LWFBrook90R_updateIsotopes_GWAT_SWAT_AdvecDiff!(u, integrator.t, integrator)
            @test p.cache_for_ADE_28[1] ≈ θ_expected
            @test u.SWATI.mm == simulation.ODESolution.u[end].SWATI.mm
            @test p.cache_for_ADE_28[14][2:end] == p.cache_for_ADE_28[18][1:end-1]
            @test p.cache_for_ADE_28[15][2:end] == p.cache_for_ADE_28[19][1:end-1]
            if uptake == :zero
                @test isnan(u.RWU.d18O) && isnan(u.RWU.d2H)
                @test (u.XYLEM.d18O, u.XYLEM.d2H) == xylem_before
            else
                @test u.RWU.d18O ≈ expected_rwu[1]
                @test u.RWU.d2H ≈ expected_rwu[2]
            end
            @test isotope_callback_allocated(integrator) == 0
        end
    end
end

@testset "Adaptive norm: unchanged fields, no temporary arrays" begin
    for isotopes in (false, true), layers in (1, 25)
        u = LWFBrook90.define_LWFB90_u0(; simulate_isotopes=isotopes,
            compute_intermediate_quantities=true, NLAYER=layers)
        u .= range(-2.0, 3.0; length=length(u))
        reference = vcat(u.GWAT.mm, u.INTS.mm, u.INTR.mm, u.SNOW.mm,
            u.CC.MJm2, u.SNOWLQ.mm, u.SWATI.mm, u.RWU.mmday, u.XYLEM.mm,
            u.TRANI.mmday, u.aux.θ, u.accum)
        expected = sqrt(sum(abs2, reference) / length(reference))
        @test LWFBrook90.norm_to_use(u, 0.0) ≈ expected
        @test norm_allocated(u) == 0
        u.aux.ψ .= NaN
        u.aux.K .= NaN
        if isotopes
            u.SWATI.d18O .= NaN
            u.GWAT.d2H = NaN
        end
        @test LWFBrook90.norm_to_use(u, 0.0) ≈ expected
        u.SWATI.mm[1] = Inf
        @test isinf(LWFBrook90.norm_to_use(u, 0.0))
        u.SWATI.mm[1] = NaN
        @test isnan(LWFBrook90.norm_to_use(u, 0.0))
        fill!(u, 0.0)
        @test LWFBrook90.norm_to_use(u, 0.0) == 0.0
    end
end
