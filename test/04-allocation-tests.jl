# Small hot-path regression checks; also runnable independently of plotting tests.
function norm_allocated(u)
    LWFBrook90.norm_to_use(u, 0.0) # compile before measuring
    return @allocated LWFBrook90.norm_to_use(u, 0.0)
end

@testset "Daily snapshots: independent accumulator values, not entire states" begin
    callbacks, _ = LWFBrook90.define_LWFB90_cb((0.0, 3.0))
    saving_callback = only(filter(cb -> hasproperty(cb.affect!, :save_func),
        callbacks.discrete_callbacks))
    save_func = saving_callback.affect!.save_func
    for isotopes in (false, true), layers in (1, 25)
        u = LWFBrook90.define_LWFB90_u0(; simulate_isotopes=isotopes,
            compute_intermediate_quantities=true, NLAYER=layers)
        u .= range(1.0, 2.0; length=length(u))
        snapshot = save_func(u, 0.0, nothing)
        expected_accum = deepcopy(u.accum)
        expected_trani = copy(u.TRANI.mmday)
        @test snapshot.accum == expected_accum
        @test propertynames(snapshot.accum) == propertynames(expected_accum)
        @test snapshot.accum.cum_d_prec == expected_accum.cum_d_prec
        @test snapshot.TRANI == expected_trani
        @test LWFBrook90.getdata(snapshot.accum) isa Vector
        @test length(LWFBrook90.getdata(snapshot.accum)) == length(u.accum)
        # The next daily reset/update must not alter previously saved fluxes.
        fill!(u.accum, 0.0)
        fill!(u.TRANI.mmday, 0.0)
        @test snapshot.accum == expected_accum
        @test snapshot.TRANI == expected_trani
        snapshot.accum.cum_d_prec = -1.0
        @test u.accum.cum_d_prec == 0.0
    end
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

# Original array-based hourly integration, retained as a numerical reference.
function inter24_reference(rain, pint, lai, sai, frintl, frints, cintrl, cintrs,
                           durations, storage, month)
    ihd = Int(floor((durations[month] + 0.1) / 2.0))
    capacity = cintrl * lai + cintrs * sai
    hourly_catch = Float64[hour < 12 - ihd || hour >= 12 + ihd ? 0.0 :
        min(1, frintl * lai + frints * sai) * rain / (2 * ihd) for hour in 0:23]
    intr = zeros(24)
    rint = zeros(24)
    irvp = zeros(24)
    intr[1] = storage
    for i in 1:24
        newint = intr[i] + (hourly_catch[i] - pint / 24)
        if newint > 0.0001
            irvp[i] = pint / 24
            rint[i] = newint > capacity ? irvp[i] + (capacity - intr[i]) : hourly_catch[i]
        else
            irvp[i] = intr[i] + hourly_catch[i]
            rint[i] = hourly_catch[i]
        end
        if i < 24
            intr[i+1] = intr[i] + (rint[i] - irvp[i])
        end
    end
    return (sum(rint), sum(irvp))
end

function plntres_reference(nlayer, soil, rtlen, relden, rtrad, rplant, fxylem, pi, rhowg)
    d = Float64.(soil.p_THICK .* (1 .- soil.p_STONEF))
    rtfrac = relden[1:nlayer] .* d ./ sum(relden[1:nlayer] .* d)
    rxylem = fxylem * rplant
    rrooti = zeros(nlayer)
    alpha = zeros(nlayer)
    for i in 1:nlayer
        if relden[i] < 0.00001 || rtlen < 0.1
            rrooti[i] = alpha[i] = 1E20
        else
            rrooti[i] = (rplant - rxylem) / rtfrac[i]
            rtdeni = rtfrac[i] * 0.001 * rtlen / d[i]
            delt = pi * rtrad^2 * rtdeni
            alpha[i] = (1 / (8 * pi * rtdeni)) * (delt - 3 - 2 * log(delt) / (1 - delt))
            alpha[i] = alpha[i] * 0.001 * rhowg / d[i]
        end
    end
    return (rxylem, rrooti, alpha)
end

function routine_allocated(f::F, args::A) where {F, A}
    f(args...) # compile before measuring
    return @allocated f(args...)
end

@testset "Daily interception: scalar hourly storage, unchanged water balance" begin
    inter24 = LWFBrook90.EVP.INTER24
    for rain in (0.0, 0.2, 4.0, 50.0), pint in (-0.2, 0.0, 1.0, 12.0),
        lai in (0.0, 4.0), storage in (0.0, 0.0001, 0.6, 5.0),
        duration in (0.0, 1.8, 4.0, 7.5, 24.0), month in (1, 12)
        # Includes dry canopy, condensation, over-capacity storage, and zero-hour storms.
        args = (rain, pint, lai, 0.5, 0.06, 0.06, 0.2, 0.1,
            fill(duration, 12), storage, month)
        actual = inter24(args...)
        expected = inter24_reference(args...)
        @test isequal(actual, expected)
        args32 = (Float32.(args[1:8])..., Float32.(args[9]), Float32(storage), month)
        @test isequal(inter24(args32...), inter24_reference(args32...))
    end
    # Two hourly-rate arrays only; allow a small margin for the returned scalar tuple.
    array_budget = routine_allocated(n -> (zeros(n), zeros(n)), (24,))
    @test routine_allocated(inter24,
        (4.0, 1.0, 4.0, 0.5, 0.06, 0.06, 0.2, 0.1, fill(4.0, 12), 0.6, 1)) <= array_budget + 64
end

@testset "Plant resistance: reuse output for root weights, preserve rootless layers" begin
    plntres = LWFBrook90.EVP.PLNTRES
    soil = (; p_THICK=[50.0, 100.0, 200.0], p_STONEF=[0.0, 0.2, 0.5])
    for relden in ([1.0, 0.5, 0.1], [1.0, 0.0, 1e-7], zeros(3)),
        rtlen in (0.0, 0.099, 0.1, 1000.0)
        args = (3, soil, rtlen, relden, 0.35, 0.125, 0.5,
            LWFBrook90.CONSTANTS.p_PI, LWFBrook90.CONSTANTS.p_RHOWG)
        roots_before = copy(relden)
        actual = plntres(args...)
        expected = plntres_reference(args...)
        @test all(a ≈ b for (a, b) in zip(actual, expected))
        @test relden == roots_before
        @test actual[2] !== relden && actual[3] !== relden
        @test routine_allocated(plntres, args) < routine_allocated(plntres_reference, args)
        soil32 = (; p_THICK=Float32.(soil.p_THICK), p_STONEF=Float32.(soil.p_STONEF))
        args32 = (3, soil32, Float32(rtlen), Float32.(relden), Float32.(args[5:end])...)
        @test isequal(plntres(args32...), plntres_reference(args32...))
    end
end
