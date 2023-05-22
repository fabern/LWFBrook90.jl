using CSV: read, File
using DataFrames: DataFrame, rename, sort!, transform!# ,select
using DataFramesMeta
using Dates: DateTime, Millisecond, Second, Day, Month, month, value, dayofyear, format
using Statistics: mean
# using Printf: @sprintf
"""
    loadSPAC(folder::String, prefix::String;
            # OPTIONAL ARGUMENTS ARE:
            simulate_isotopes::Bool = true,
            canopy_evolution  = "meteoveg.csv",
            Δz_thickness_m    = "soil_discretization.csv",
            root_distribution = "soil_discretization.csv",
            IC_soil           = "soil_discretization.csv",
            IC_scalar         = "initial_conditions.csv",
            storm_durations_h = "meteo_storm_durations.csv")

Create instance of SPAC model by loading different input files for LWFBrook90 from folder `folder`` with the names:
- `[PREFIX]_meteoveg.csv`
- `[PREFIX]_meteo_iso.csv` (needed for isotopic)
- `[PREFIX]_param.csv`
- `[PREFIX]_meteo_storm_durations.csv`
- `[PREFIX]_initial_conditions.csv`
- `[PREFIX]_soil_horizons.csv`

These files were created with an R script `generate_LWFBrook90jl_Input.R` that
takes the same arguements as the R function `LWFBrook90R::run_LWFB90()` and generates
the corresponding input files for LWFBrook90.jl.

Soil discretization can be provided as vector with the thickness of each cell, e.g: `Δz_thickness_m = fill(0.04, 20)`.

Root distribution can be provided as arguments to function `Rootden_beta_()` as  `(beta = 0.98, )` or `(beta = 0.98, z_rootMax_m=-0.5)`.

Meteo storm durations can be provided as vector for each month, e.g.: `[4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]`

Initial conditions in the soil can be provided as NamedTuple `IC_soil = (PSIM_init_kPa = -7.0, delta18O_init_permil = -9.0, delta2H_init_permil = -11.0)`

Canopy evolution can be provided as NamedTuple, giving the constant values of DENSEF, HEIGHT, SAI and dynamic evolution of LAI:

    canopy_evolution = (DENSEF_rel = 100,
                        HEIGHT_rel = 100,
                        SAI_rel    = 100,
                        LAI_rel = (DOY_Bstart = 120,
                                Bduration  = 21,
                                DOY_Cstart = 270,
                                Cduration  = 60,
                                LAI_perc_BtoC = 100,
                                LAI_perc_CtoB = 20))

The vegetation season in terms of LAI is defined by defining the year into 4 parts:
A->B, B->B+, B+->C, C->C+, C+->A where A is the start of the year (January 1st).
effectively giving the 4 parts: C+->B, B->B+, B+->C, C->C+.
The position and durationof these periods are defined in days by parameters: DOY\\_Bstart, Bduration, DOY\\_Cstart, and Cduration.
LAI (in percent) is constant from C+->B and B+->C (given by LAI\\_perc\\_CtoB and LAI\\_perc\\_BtoC)
and a simple linear interpolation is done for in\\_between (i.e. budburst and leaffall).

Initial conditions of states other than soil can be provided as NamedTuple, e.g.:

    IC_scalar = (amount = (u_GWAT_init_mm = 0,
                           u_INTS_init_mm = 0,
                           u_INTR_init_mm = 0,
                           u_SNOW_init_mm = 0,
                           u_CC_init_MJ_per_m2 = 0,
                           u_SNOWLQ_init_mm =  0),
                d18O    = (u_GWAT_init_permil = -13.,
                           u_INTS_init_permil = -13.,
                           u_INTR_init_permil = -13.,
                           u_SNOW_init_permil = -13.),
                d2H     = (u_GWAT_init_permil = -95.,
                           u_INTS_init_permil = -95.,
                           u_INTR_init_permil = -95.,
                           u_SNOW_init_permil = -95.))
"""
function loadSPAC(folder::String, prefix::String;
    simulate_isotopes::Bool = true,
    # compute_intermediate_quantities::Bool = true,
    canopy_evolution  = "meteoveg.csv",
    Δz_thickness_m    = "soil_discretization.csv",
    root_distribution = "soil_discretization.csv",  #either "soil_discretization.csv" or (beta = 0.97, z_rootMax_m = nothing)
    IC_soil           = "soil_discretization.csv",  #either "soil_discretization.csv" or (PSIM_init_kPa = -6.3, delta18O_init_permil = -10., delta2H_init_permil = -95.)
    IC_scalar         = "initial_conditions.csv",
        #either "initial_conditions.csv" or IC_scalar = (amount = (u_GWAT_init_mm = 0, ...),
        #                                                d18O   = (u_GWAT_init_permil = -13., ...),
        #                                                d2H    = (u_GWAT_init_permil = -95., ...))
    storm_durations_h = "meteo_storm_durations.csv")  #either "meteo_storm_durations.csv" or [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]

    ## Define paths of all input files
    input_file_XXXX = prefix*"_XXXX"*".csv"
    path_meteoveg           = joinpath(folder, replace(input_file_XXXX, "XXXX" => "meteoveg"))
    path_param              = joinpath(folder, replace(input_file_XXXX, "XXXX" => "param"))
    # unused path_precdat   = joinpath(folder, replace(input_file_XXXX, "XXXX" => "precdat"))
    path_storm_durations    = joinpath(folder, replace(input_file_XXXX, "XXXX" => "meteo_storm_durations"))
    path_initial_conditions = joinpath(folder, replace(input_file_XXXX, "XXXX" => "initial_conditions"))
    path_soil_horizons      = joinpath(folder, replace(input_file_XXXX, "XXXX" => "soil_horizons"))
    path_soil_discretization = joinpath(folder, replace(input_file_XXXX, "XXXX" => "soil_discretization"))

    ## Define solver options
    solver_options =
        (#Reset                           = false, # currently only Reset = 0 implemented
         compute_intermediate_quantities = true,   # Flag whether ODE containes additional quantities than only states
         simulate_isotopes               = simulate_isotopes
        )

    ## Load model input parameters
    params, solver_opts = init_param(path_param; simulate_isotopes = simulate_isotopes)

    solver_options = merge(solver_options, solver_opts) # append to manually provided solver options

    ## Load time-varying atmospheric forcing
    reference_date, input_meteoveg, meteo_iso_forcing, storm_durations =
        init_forcing(path_meteoveg, path_storm_durations; simulate_isotopes, storm_durations_h)

    tspan = extrema(input_meteoveg.days)

    meteo_forcing = input_meteoveg[:, [:days, :GLOBRAD, :TMAX, :TMIN, :VAPPRES, :WIND, :PRECIN]]
    @assert all(meteo_forcing.GLOBRAD .>= 0) "Error in vegetation parameters: GLOBRAD must be above 0."
    # @assert all(meteo_forcing.TMAX    .> 0) "Error in vegetation parameters: TMAX must be above 0."
    # @assert all(meteo_forcing.TMIN    .> 0) "Error in vegetation parameters: TMIN must be above 0."
    @assert all(meteo_forcing.VAPPRES .>= 0) "Error in vegetation parameters: VAPPRES must be above 0."
    @assert all(meteo_forcing.WIND    .>= 0) "Error in vegetation parameters: WIND must be above 0."
    @assert all(meteo_forcing.PRECIN  .>= 0) "Error in vegetation parameters: PRECIN must be above 0."

    ## Load time-varying vegetation parameters
    if (canopy_evolution == "meteoveg.csv")
        # Use DataFrame from meteoveg.csv
        if ("DENSEF_rel" in names(input_meteoveg) && "HEIGHT_rel" in names(input_meteoveg) &&
            "LAI_rel"    in names(input_meteoveg) && "SAI_rel"    in names(input_meteoveg))
            _to_use_canopy_evolution = input_meteoveg[:, [:days, :DENSEF_rel, :HEIGHT_rel, :LAI_rel, :SAI_rel]] # store DataFrame in SPAC
            # Assert validity of vegetation values
            @assert all(_to_use_canopy_evolution.LAI_rel .>= 0) "Error in vegetation parameters: LAI must be above 0%."
        else
            error("""
                Input_meteoveg is expected to contain one or multiple of the columns: :DENSEF_re, :HEIGHT_rel, :LAI_rel, or :SAI_rel.
                Please check your input files with the current documentation and possibly contact the developer team if the error persists.
                If it does not contain the columns please provide a parametrization to loadSPAC(, canopy_evolution::NamedTuple = ())
                with the NamedTuple containing relative values in percent:
                `(DENSEF_rel = 100, HEIGHT_rel = 100, SAI_rel = 100, LAI_rel = (DOY_Bstart, Bduration, DOY_Cstart, Cduration, LAI_perc_BtoC, LAI_perc_CtoB))`
                """)
        end
    else
        # Use received parameter from arguments to loadSPAC()
        if ("DENSEF_rel" in names(input_meteoveg) && "HEIGHT_rel" in names(input_meteoveg) &&
            "LAI_rel"    in names(input_meteoveg) && "SAI_rel"    in names(input_meteoveg))
            @warn "Received canopy_evolution in loadSPAC(), overwriting values from `meteoveg.csv`."
        end
        @assert canopy_evolution isa NamedTuple
        @assert keys(canopy_evolution) == (:DENSEF_rel, :HEIGHT_rel, :SAI_rel, :LAI_rel)
        @assert keys(canopy_evolution.LAI_rel) == (:DOY_Bstart, :Bduration, :DOY_Cstart, :Cduration, :LAI_perc_BtoC, :LAI_perc_CtoB)

        _to_use_canopy_evolution = canopy_evolution # store parameter arguments in SPAC
        # Assert validity of vegetation values
        @assert all(_to_use_canopy_evolution.LAI_rel.LAI_perc_BtoC >= 0) "Error in vegetation parameters: LAI must be above 0%."
        @assert all(_to_use_canopy_evolution.LAI_rel.LAI_perc_CtoB >= 0) "Error in vegetation parameters: LAI must be above 0%."

        # NOTE: that values in canopy_evolution have to be relative express in percent
        #       Absolute values are derived in combination with:
        #           params.DENSEF_baseline_
        #           params.SAI_baseline_
        #           params.AGE_baseline_yrs
        #           params.HEIGHT_baseline_m
        #           params.MAXLAI
    end

    # Assert validity of vegetation values
    @assert all(_to_use_canopy_evolution.DENSEF_rel .> 5) "DENSEF (in meteoveg.csv or canopy_evolution argument) should not be set lower than 5% as it affects aerodynamics."
    @assert all(_to_use_canopy_evolution.HEIGHT_rel .> 0) "Error in vegetation parameters: HEIGHT must be above 0%."
    @assert all(_to_use_canopy_evolution.SAI_rel .> 0) "Error in vegetation parameters: SAI must be above 0%."

    ## Load space-varying soil data
    soil_horizons = init_soil(path_soil_horizons)

    ## Load initial conditions of scalar state variables
    if IC_scalar isa NamedTuple
        if isfile(path_initial_conditions) @warn "Requested to overwrite initial conditions. Values in $path_initial_conditions are ignored." end
        IC_scalar =
            DataFrame(u_GWAT_init_mm      = [IC_scalar.amount.u_GWAT_init_mm     ,IC_scalar.d18O.u_GWAT_init_permil ,IC_scalar.d2H.u_GWAT_init_permil],
                      u_INTS_init_mm      = [IC_scalar.amount.u_INTS_init_mm     ,IC_scalar.d18O.u_INTS_init_permil ,IC_scalar.d2H.u_INTS_init_permil],
                      u_INTR_init_mm      = [IC_scalar.amount.u_INTR_init_mm     ,IC_scalar.d18O.u_INTR_init_permil ,IC_scalar.d2H.u_INTR_init_permil],
                      u_SNOW_init_mm      = [IC_scalar.amount.u_SNOW_init_mm     ,IC_scalar.d18O.u_SNOW_init_permil ,IC_scalar.d2H.u_SNOW_init_permil],
                      u_CC_init_MJ_per_m2 = [IC_scalar.amount.u_CC_init_MJ_per_m2,NaN                               ,NaN                             ],
                      u_SNOWLQ_init_mm    = [IC_scalar.amount.u_SNOWLQ_init_mm   ,NaN                               ,NaN                             ])
    elseif IC_scalar == "initial_conditions.csv"
        IC_scalar = init_IC(path_initial_conditions)
    else
        error("Unknown format for argument `IC_scalar`: $IC_scalar")
    end

    ## Load soil discretizations either from `soil_discretizations.csv`
    ## or then use provided Δz_thickness_m which also requires IC and root distribution parameters

    # Possible cases of arguments:
    # (Δz_thickness_m="soil_discretization.csv", root_distribution="soil_discretization.csv", IC_soil="soil_discretization.csv"): load soil_discretization.csv
    # (Δz_thickness_m="soil_discretization.csv", root_distribution="soil_discretization.csv", IC_soil=(.=, .=, .=)             ): load soil_discretization.csv + overwrite IC
    # (Δz_thickness_m="soil_discretization.csv", root_distribution=(.=, .=, .=),              IC_soil="soil_discretization.csv"): load soil_discretization.csv                + overwrite roots
    # (Δz_thickness_m="soil_discretization.csv", root_distribution=(.=, .=, .=),              IC_soil=(.=, .=, .=)             ): load soil_discretization.csv + overwrite IC + overwrite roots

    # (Δz_thickness_m=[., .]                   , root_distribution="soil_discretization.csv", IC_soil="soil_discretization.csv"): # error: provide parametric version of root distribution and IC
    # (Δz_thickness_m=[., .]                   , root_distribution="soil_discretization.csv", IC_soil=(.=, .=, .=)             ): # error: provide parametric version of root distribution
    # (Δz_thickness_m=[., .]                   , root_distribution=(.=, .=, .=),              IC_soil="soil_discretization.csv"): # error: provide parametric version of                       IC
    # (Δz_thickness_m=[., .]                   , root_distribution=(.=, .=, .=),              IC_soil=(.=, .=, .=)             ): # okay. warning: not using "soil_discretization.csv"

    if Δz_thickness_m isa Vector
        # Define soil discretiztion freely
        if isfile(path_soil_discretization) @warn "loadSPAC(...; Δz_thickness_m = ...) provided. Ignoring `soil_discretization.csv`, i.e. also overwriting root distribution and initial conditions." end
        interfaces_m = -vcat(0, cumsum(Δz_thickness_m))
        soil_discretization = DataFrame(Upper_m           = interfaces_m[Not(end)],
                                    Lower_m               = interfaces_m[Not(1)],
                                    Rootden_              = NaN,
                                    uAux_PSIM_init_kPa    = NaN,
                                    u_delta18O_init_permil= NaN,
                                    u_delta2H_init_permil = NaN)

        if root_distribution == "soil_discretization.csv" error("Requested to create soil discretization manually instead of using `soil_discretization.csv`, but no parametric root_distribution provided.") end
        if IC_soil           == "soil_discretization.csv" error("Requested to create soil discretization manually instead of using `soil_discretization.csv`, but no parametric soil initial conditions provided.") end
        # _to_use_root_distribution = (beta = 0.97, z_rootMax_m = nothing)
        # _to_use_IC_soil           = (PSIM_init_kPa = -6.3, delta18O_init_permil = -10., delta2H_init_permil = -95.)
        _to_use_Δz_thickness_m = Δz_thickness_m

        # Overwrite soil_discretization with roots # TODO: make this a function to reuse if needed in setup()
        _to_use_root_distribution = root_distribution
        overwrite_rootden!(soil_discretization, _to_use_root_distribution, _to_use_Δz_thickness_m)
        # Overwrite soil_discretization with IC # TODO: make this a function to reuse if needed in setup()
        _to_use_IC_soil = IC_soil
        overwrite_IC!(soil_discretization, _to_use_IC_soil, simulate_isotopes)

    elseif Δz_thickness_m == "soil_discretization.csv"
        if isfile(path_soil_discretization)
            soil_discretization = LWFBrook90.read_path_soil_discretization(path_soil_discretization)
            # Assert type stability by executing disallowmissing!
            # Impose type of Float64 instead of Float64?, by defining unused variables as -9999.99
            if !("u_delta18O_init_permil" in names(soil_discretization)) insertcols!(soil_discretization, :u_delta18O_init_permil => -9999.0) end #NaN) end
            if !("u_delta2H_init_permil"  in names(soil_discretization)) insertcols!(soil_discretization, :u_delta2H_init_permil => -9999.0)  end #NaN) end
            if (any(ismissing.(soil_discretization.u_delta18O_init_permil)) || any(ismissing.(soil_discretization.u_delta2H_init_permil)))
                soil_discretization.u_delta18O_init_permil .= -9999.0
                soil_discretization.u_delta2H_init_permil  .= -9999.0
            end
            disallowmissing!(soil_discretization, [:Rootden_, :uAux_PSIM_init_kPa, :u_delta18O_init_permil, :u_delta2H_init_permil])

            _to_use_Δz_thickness_m = soil_discretization.Upper_m - soil_discretization.Lower_m
            _to_use_root_distribution = root_distribution
            _to_use_IC_soil = IC_soil

            if root_distribution != "soil_discretization.csv"
                # overwrite root distribution and warn about overwriting
                overwrite_rootden!(soil_discretization, _to_use_root_distribution, _to_use_Δz_thickness_m)
                @warn "Requested to overwrite root distribution. Root distribution defined in soil_discretization is ignored."
            end
            if IC_soil != "soil_discretization.csv"
                # overwrite initial conditions and warn about overwriting
                overwrite_IC!(soil_discretization, _to_use_IC_soil, simulate_isotopes)
                @warn "Requested to overwrite initial conditions. Initial conditions defined in soil_discretization are ignored."
            end
        else
            error("No file `$path_soil_discretization` found. Either define a file `soil_discretiztions.csv` or then provide loadSPAC(...; soil_discretization = [0.04, 0.04, 0.04, 0.04, 0.04].")
            # @warn """
            # No file `$path_soil_discretization` found.
            # soil_discretization is derived from layers in `$path_soil_horizons`.
            # Consider providing either a file `soil_discretiztions.csv` or then provide loadSPAC(...; soil_discretization = [0.04, 0.04, 0.04, 0.04, 0.04])
            # Using default initial conditions (-6.3 kPa, ... permil, ... permil) and root distribution (beta = 0.97), unless overwritten in setup().
            # """
            # soil_discretization = DataFrame(Upper_m                = soil_horizons.Upper_m,
            #                             Lower_m               = soil_horizons.Lower_m,
            #                             Rootden_              = NaN,
            #                             uAux_PSIM_init_kPa    = NaN,
            #                             u_delta18O_init_permil= NaN,
            #                             u_delta2H_init_permil = NaN)
            # # default root distribution and initial conditiosn unless otherwise defined in setup(), then they are overwritten
            # root_distribution = (beta = 0.97, z_rootMax_m = nothing)
            # IC_soil           = (PSIM_init_kPa = -6.3, delta18O_init_permil = -10., delta2H_init_permil = -95.)
        end
    else
        error("Unknown format for argument `Δz_thickness_m`: $Δz_thickness_m")
    end

    ## Extend soil horizons if needed by requested soil discretization
    # (in such a case emit a warning)
    extended_soil_horizons = extend_lowest_horizon(soil_horizons, soil_discretization)

    ## Make time dependent input parameters continuous in time (interpolate)
    (meteo_forcing_cont, meteo_iso_forcing_cont) =
        LWFBrook90.interpolate_meteo(;
            meteo_forcing                 = meteo_forcing,
            meteo_iso_forcing             = meteo_iso_forcing);

    return SPAC(;
        reference_date    = reference_date,
        tspan             = tspan,
        solver_options    = solver_options,
        soil_discretization    = (Δz = _to_use_Δz_thickness_m,
                                  df = soil_discretization),
        forcing = (meteo           = meteo_forcing_cont,
                   meteo_iso       = meteo_iso_forcing_cont,
                   storm_durations = storm_durations),
        pars    = (params = params,
                   root_distribution = _to_use_root_distribution,
                   IC_scalar = IC_scalar,
                   IC_soil   = _to_use_IC_soil,
                   canopy_evolution = _to_use_canopy_evolution,
                   soil_horizons = extended_soil_horizons),
        )
end

function Base.show(io::IO, ::MIME"text/plain", models::Vector{DiscretizedSPAC}) print(io, "Array of $(length(models)) DiscretizedSPAC. Not showing details...") end # https://discourse.julialang.org/t/9589
function Base.show(io::IO, mime::MIME"text/plain", model::SPAC; show_SPAC_title=true)
    if (show_SPAC_title) println(io, "SPAC model:") end

    println(io, "===== DATES:===============")
    println(io, format(model.reference_date, "YYYY-mm-dd"), ", tspan:", model.tspan, "days")
    println(io, "\n===== METEO FORCING:===============")
    show_avg_and_range = function(vector, title)
        # "$(title)avg:$(round(Statistics.mean(vector); digits=2)), range:$(extrema(vector))"
        @sprintf("%savg:%7.2f, range:%7.2f to%7.2f",
                 title, mean(vector), extrema(vector)[1], extrema(vector)[2])
    end
    println(io, show_avg_and_range(model.forcing.meteo["p_GLOBRAD"].itp.itp.coefs,     "GLOBRAD (MJ/m2/day): "))
    println(io, show_avg_and_range(model.forcing.meteo["p_PREC"].itp.itp.coefs,        "PREC       (mm/day): "))
    println(io, show_avg_and_range(model.forcing.meteo["p_TMAX"].itp.itp.coefs,        "TMAX           (°C): "))
    println(io, show_avg_and_range(model.forcing.meteo["p_TMIN"].itp.itp.coefs,        "TMIN           (°C): "))
    println(io, show_avg_and_range(model.forcing.meteo["p_VAPPRES"].itp.itp.coefs,     "VAPPRES       (kPa): "))
    println(io, show_avg_and_range(model.forcing.meteo["p_WIND"].itp.itp.coefs,        "WIND          (m/s): "))
    if (model.solver_options.simulate_isotopes)
        println(io, show_avg_and_range(model.forcing.meteo_iso["p_d18OPREC"].itp.coefs,"δ18O            (‰): "))
        println(io, show_avg_and_range(model.forcing.meteo_iso["p_d2HPREC"].itp.coefs, "δ2H             (‰): "))
    end
    # println(io, model.forcing.storm_durations)

    println(io, "\n===== CANOPY EVOLUTION:===============")
    if model.pars.canopy_evolution isa DataFrame
        println(io, "model.pars.canopy_evolution was loaded form meteoveg.csv")
        # println(io, show_avg_and_range(model.pars.canopy_evolution.DENSEF_rel,"DENSEF            (%): "))
        # println(io, show_avg_and_range(model.pars.canopy_evolution.HEIGHT_rel,"HEIGHT            (%): "))
        # println(io, show_avg_and_range(model.pars.canopy_evolution.LAI_rel,"LAI            (%): "))
        # println(io, show_avg_and_range(model.pars.canopy_evolution.SAI_rel,"SAI            (%): "))
    else
        println(io, model.pars.canopy_evolution)
    end

    println(io, "\n===== INITIAL CONDITIONS:===============")
    println(io, "Soil   IC: $(model.pars.IC_soil)")
    print(  io, "Scalar IC: ")
    println(io, model.pars.IC_scalar)

    println(io, "\n===== MODEL PARAMETRIZATION:===============")
    # println(io, model.params)
    #display(io, model.params) # dump(model.params); using PrettyPrinting; pprintln(model.params)
    # maxlength = maximum(length.(string.(keys(model.params)))) # 17
    # maxlength # 17
    # for (k,v) in pairs(model.params) println(@sprintf("%18s",k) => v) end
    # for (k,v) in pairs(model.params) println(@sprintf("%18s",k) => @sprintf("% 8.1f",v)) end
    # for (k,v) in pairs(model.params) println(@sprintf("%18s => % 8.1f",k,v)) end
    string_vec = [@sprintf("%18s => % 8.1f",k,v) for (k,v) in pairs(model.pars.params)];
    # reshape(string_vec, 27, 3)
    # reduce.(*, reshape(string_vec, 27, 3))
    # reduce.(*, reshape(vcat(string_vec, fill("", 21*4-length(string_vec))), 21, 4))
    # reshape(vcat(string_vec, fill("", 21*4-length(string_vec))), 21, 4)
    show(IOContext(io, :limit => true), "text/plain",
         reshape(vcat(string_vec, fill("", 21*4-length(string_vec))), 21, 4))
    # show(IOContext(io, :limit => true), "text/plain",
    #      reshape(vcat(string_vec, fill("", 42*2-length(string_vec))), 42, 2))

    println(io, "\n===== SOIL DOMAIN:===============")
    print(  io, "Root distribution:       "); println(io, model.pars.root_distribution)
    print(  io, "Soil layer properties:   ")
    # println(io, model.pars.soil_horizons)
    # for shp in model.pars.soil_horizons.shp println(io, shp) end
    # show(IOContext(io, :limit => false), "text/plain", model.pars.soil_horizons)
    show(IOContext(io, :limit => false), mime, model.pars.soil_horizons[:,Not(:HorizonNr)]; truncate = 100);
    println(io, "")
    Δz = model.soil_discretization.Δz
    println(io, "Soil discretized into N=$(length(model.soil_discretization.df.Lower_m)) layers, "*
                "$(@sprintf("Δz layers: (avg, min, max) = (%.3f,%.3f,%.3f)m.", mean(Δz),minimum(Δz),maximum(Δz)))")
    println(io, round.(model.soil_discretization.df.Lower_m; digits=3))

    println(io, "\n===== SOLVER OPTIONS:===============")
    println(io, model.solver_options)
end

######################
# Define functions to handle DateTimes and convert into Days as Floats
"""
    DateTime2RelativeDaysFloat(x,reference_DateTime)

Transforms DateTimes `x` to simulation time in days since reference date.
"""
function DateTime2RelativeDaysFloat(x::DateTime, reference::DateTime)
    ms2days = 1.0/(24.0*3600.0*1000.0) # to convert milliseconds to days
    ms2days*value(convert(Millisecond, x-reference))
end
"""
    RelativeDaysFloat2DateTime(t, reference_DateTime)

Transforms simulation time `t` in days since reference date to DateTimes.
"""
function RelativeDaysFloat2DateTime(t::Float64, reference::DateTime)
    # reference + Day(floor(t))
    t_sec = 60*60*24*t # t is in days, t_sec in seconds
    reference + Second(floor(t_sec))
end
"""
    p_DOY(t::Float64, reference::DateTime)

Get DOY (Day Of Year) from simulation time
"""
function p_DOY(t::Float64, reference::DateTime)
    dayofyear(reference + Day(floor(t)))
end
"""
    p_MONTHN(t::Float64, reference::DateTime)

Get Month from simulation time
"""
function p_MONTHN(t::Float64, reference::DateTime)
    month(reference + Day(floor(t)))
end

function init_forcing(path_meteoveg, path_storm_durations; simulate_isotopes = true, storm_durations_h)

    # Load daily values of meteo and precipitation isotopes
    if (simulate_isotopes)
        path_meteoiso = replace(path_meteoveg, "meteoveg" => "meteoiso")
    end
    meteo_forcing, reference_date = read_path_meteoveg(path_meteoveg)

    if (simulate_isotopes)
        meteo_iso_forcing, meteo_forcing = read_path_meteoiso(path_meteoiso,
            meteo_forcing,
            reference_date)
    else
        meteo_iso_forcing = nothing
    end

    # Load parameter for subdaily model of precipitation
    #' @param storm_durations a [1,12]-matrix of precipitation durations (hours) for each month.
    if storm_durations_h isa Vector
        @assert length(storm_durations_h) == 12 "Wrong format for argument `storm_durations_h`: $storm_durations_h. It must be a vector of length 12."
        storm_durations = DataFrame(
            month             = ["January", "Februrary", "March", "April", "May", "June",
                                 "July", "August", "September", "October", "November", "December"],
            storm_durations_h = storm_durations_h)
            # storm_durations_h = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0])
        if isfile(path_storm_durations) @warn "Requested to overwrite storm durations. Values in $path_storm_durations are ignored." end
    elseif storm_durations_h == "meteo_storm_durations.csv"
        storm_durations = read_path_storm_durations(path_storm_durations)
    else
        error("Unknown format for argument `meteo_storm_durations`: $storm_durations_h")
    end

    # Impose type of Float64 instead of Float64?
    disallowmissing!(meteo_forcing)
    if (simulate_isotopes)
        disallowmissing!(meteo_iso_forcing, [:days, :delta18O_permil, :delta2H_permil])
    end
    disallowmissing!(storm_durations, [:month, :storm_durations_h])

    return reference_date, meteo_forcing, meteo_iso_forcing, storm_durations
end

function init_IC(path_initial_conditions; simulate_isotopes = true)
    read_path_initial_conditions(path_initial_conditions)
end

function init_soil(path_soil_horizons)
    read_path_soil_horizons(path_soil_horizons)
end

function init_param(path_param; simulate_isotopes = true)
    input_param, solver_opts = read_path_param(path_param; simulate_isotopes = simulate_isotopes)
end

# path_meteoveg = "examples/BEA2016-reset-FALSE-input/BEA2016-reset-FALSE_meteoveg.csv"
function read_path_meteoveg(path_meteoveg)

    if (File(path_meteoveg).cols == 7)
        # [:dates, :globrad_MJDayM2, :tmax_degC, :tmin_degC, :vappres_kPa, :windspeed_ms, :prec_mmDay]
        # accomodate case when canopy_evolution is provided manually:
        parsing_types =
            Dict(:dates          => DateTime,
                :globrad_MJDayM2 => Float64,
                :tmax_degC       => Float64,
                :tmin_degC       => Float64,
                :vappres_kPa     => Float64,
                :windspeed_ms    => Float64,
                :prec_mmDay      => Float64)
    else
        parsing_types =
            Dict(:dates          => DateTime,
                :globrad_MJDayM2 => Float64,
                :tmax_degC       => Float64,
                :tmin_degC       => Float64,
                :vappres_kPa     => Float64,
                :windspeed_ms    => Float64,
                :prec_mmDay      => Float64,
                :densef_percent  => Float64,
                :height_percent  => Float64,
                :lai_percent     => Float64,
                :sai_percent     => Float64)
    end
    # input_meteoveg = DataFrame(File(path_meteoveg;
    #     skipto=3, delim=',', ignorerepeated=false,
    #     # Be strict about loading NA's -> error if NA present
    #     types=parsing_types, missingstring = nothing, strict=true))
    # input_meteoveg.dates = DateTime.(input_meteoveg.dates)
    # ipnut_meteove = transform(input_meteove, :dates = DateTime.(:dates))
    input_meteoveg = @chain begin DataFrame(File(path_meteoveg;
        skipto=3, delim=',', ignorerepeated=false,
        # Be strict about loading NA's -> error if NA present
        types=parsing_types, missingstring = nothing, strict=true))
        transform(:dates => (d) -> DateTime.(d), renamecols = false)
    end
    expected_names = [String(k) for k in keys(parsing_types)]
    assert_colnames_as_expected(input_meteoveg, path_meteoveg, expected_names)
    # DataFrame(dates = "YYYY-MM-DD", globrad_MJDayM2 = "MJ/Day/m2",
    #     tmax_degC = "degree C", tmin_degC = "degree C", vappres_kPa = "kPa",
    #     windspeed_ms = "m per s", prec_mmDay = "mm per day",
    #     densef_percent = "percent", height_percent = "percent",
    #     lai_percent = "percent", sai_percent = "percent")

    # Assert units:
    assert_unitsHeader_as_expected(path_meteoveg,
        DataFrame(dates = "YYYY-MM-DD", globrad_MJDayM2 = "MJ/Day/m2",
        tmax_degC = "degree C", tmin_degC = "degree C", vappres_kPa = "kPa",
        windspeed_ms = "m per s", prec_mmDay = "mm per day",
        densef_percent = "percent", height_percent = "percent",
        lai_percent = "percent", sai_percent = "percent")[[1], 1:ncol(input_meteoveg)])

    # Assert validity of values
    # ...

    # Assert that no gaps
    @assert all(diff(input_meteoveg.dates) .== Millisecond(86400000)) "There are gaps in the forcing file. The file ($path_meteoveg) needs to have a value for each data from start until the end."

    # Identify period of interest
    # Starting date: latest among the input data
    # Stopping date: earliest among the input data
    starting_date = maximum(minimum,[input_meteoveg[:,"dates"]])
    stopping_date = minimum(maximum,[input_meteoveg[:,"dates"]])

    input_meteoveg = filter(:dates => (d) -> d .>= starting_date && d .<= stopping_date, input_meteoveg)

    if (File(path_meteoveg).cols == 7)
        rename!(input_meteoveg,
            :globrad_MJDayM2 => :GLOBRAD,
            :tmax_degC       => :TMAX,
            :tmin_degC       => :TMIN,
            :vappres_kPa     => :VAPPRES,
            :windspeed_ms    => :WIND,
            :prec_mmDay      => :PRECIN)
    else
        rename!(input_meteoveg,
            :globrad_MJDayM2 => :GLOBRAD,
            :tmax_degC       => :TMAX,
            :tmin_degC       => :TMIN,
            :vappres_kPa     => :VAPPRES,
            :windspeed_ms    => :WIND,
            :prec_mmDay      => :PRECIN,
            :densef_percent  => :DENSEF_rel,
            :height_percent  => :HEIGHT_rel,
            :lai_percent     => :LAI_rel,
            :sai_percent     => :SAI_rel)
    end

    # Transform times from DateTimes to simulation time (Float of Days)
    reference_date = starting_date

    input_meteoveg = @chain input_meteoveg begin
        transform(:dates => (d) -> DateTime2RelativeDaysFloat.(d, reference_date), renamecols = false)
        rename(Dict(:dates => :days))
        end

    return input_meteoveg, reference_date
end

function read_path_meteoiso(path_meteoiso,
        input_meteoveg,
        input_meteoveg_reference_date)

    #####
    # If Data contains as first line: "^#Data produced with Piso.AI"
    # then apply special treatment of input file, assuming it comes directly as a download from:
    # https://isotope.bot.unibas.ch/PisoAI/ or https://isotope.bot.unibas.ch/PisoAI-eur1900-v1-2020/

    # Special treatment of Piso.AI data:
    # Piso.AI (or GNIP) attributes data to the 15th of each month but is
    # representative for the monthly data.
    # In that case: shift dates to the end of the month (and repeat the very first
    # measurement) in order to interpolate piecewise constant backward in time.

    # Otherwise assume dates reported are the end dates of the collection period.
    # In that case, extend the first measurement period to the beginning of the simulation
    # and interpolate constant backward in time.

    is_from_PisoAI = occursin(r"^#Data produced with Piso.AI",readlines(path_meteoiso)[1])

    if is_from_PisoAI
        parsing_types =
                Dict(#:Column1        => Int64, # is directly removed afterwards
                    :Site            => Char,
                    :Date            => DateTime,
                    :Latitude        => Float64,
                    :Longitude       => Float64,
                    :Elevation       => Float64,
                    Symbol("d18O.Piso.AI")     => Float64,
                    Symbol("d2H.Piso.AI")      => Float64)

        input_meteoiso = @chain begin DataFrame(File(path_meteoiso; header = 4,
                    skipto=5, delim=',', ignorerepeated=true,
                    # Don't be strict, allow for NA as missing. Treat this later with disallowmissing!.
                    types=parsing_types, missingstring = ["","NA"]))
                    end
        select!(input_meteoiso, Not([:Column1]))
    else
        # else the dataframe should be of the following structure:
            # dates,delta18O_permil,delta2H_permil
            # YYYY-MM-DD,permil,permil
            # 2021-01-01,NA,NA
            # 2021-01-15,-10.1,-70.1
            # 2021-01-29,-10.1,-70.1
            # ...
        # where the first line
        parsing_types =
                Dict(:dates          => DateTime,
                    :delta18O_permil => Float64,
                    :delta2H_permil  => Float64)

        input_meteoiso = DataFrame(File(path_meteoiso;
                skipto=3, delim=',', ignorerepeated=true,
                # Don't be strict, allow for NA as missing. Treat this later with disallowmissing!.
                types=parsing_types, missingstring = ["","NA"]))
    end

    # Assert column names as expected
    expected_names = [String(k) for k in keys(parsing_types)]
    assert_colnames_as_expected(input_meteoiso, path_meteoiso, expected_names)

    if is_from_PisoAI
        # -> a) only keep needed columns and rename them
        rename!(input_meteoiso,
                :Date                  => :dates,
                Symbol("d18O.Piso.AI") => :delta18O_permil,
                Symbol("d2H.Piso.AI")  => :delta2H_permil)
        select!(input_meteoiso, [:dates, :delta18O_permil, :delta2H_permil])
        # select!(input_meteoiso, Not([:Site, :Latitude, :Longitude, :Elevation]))

        # -> b) Shift all measured dates to end of month to be able to constant interpolate backwards.
        #       Further repeat the very line concerning the first month in the timeseries to be
        #       able to interpolate between the first and last day of the first month.
        input_meteoiso_first = deepcopy(input_meteoiso[1,:])
        input_meteoiso_first.dates = floor(input_meteoiso_first.dates, Month)
        input_meteoiso = @chain input_meteoiso begin
                transform(:dates => (d) -> ceil.(d, Month), renamecols = false)
                end
        push!(input_meteoiso, input_meteoiso_first) # repeat the first
        sort!(input_meteoiso, [:dates])

        # Modify last to remain within year (31.12. instead of 01.01)
        input_meteoiso[end,:].dates = input_meteoiso[end,:].dates - Day(1)
    else
        # Assert units as expected
        assert_unitsHeader_as_expected(path_meteoiso,
            DataFrame(dates="YYYY-MM-DD",delta18O_permil="permil",delta2H_permil="permil"))

        # Repeat the value of the second row to the first row
        # (Note, that the first row represents the start date of the period reported in the
        # second row)
        input_meteoiso[1,:delta18O_permil] = input_meteoiso[2,:delta18O_permil]
        input_meteoiso[1,:delta2H_permil]  = input_meteoiso[2,:delta2H_permil]
    end

    # Impose type of Float64 instead of Float64?
    disallowmissing!(input_meteoiso, [:dates, :delta18O_permil, :delta2H_permil])
    #####

    # Transform times from DateTimes to simulation time (Float of Days)
    input_meteoiso = @chain input_meteoiso begin
        transform(:dates => (d) -> DateTime2RelativeDaysFloat.(d, input_meteoveg_reference_date), renamecols = false)
        rename(Dict(:dates => :days))
        end

    # Check period
    # Check if overlap with meteoveg exists
    startday_iso = minimum(input_meteoiso[:,"days"])
    startday_veg = minimum(input_meteoveg[:,"days"])
    stopday_iso  = maximum(input_meteoiso[:,"days"])
    stopday_veg  = maximum(input_meteoveg[:,"days"])
    startday = max(startday_iso, startday_veg)
    endday   = min(stopday_iso, stopday_veg)

    if ((stopday_iso > stopday_veg) | (startday_iso < startday_veg))
        @warn """
        Isotopic signature of precipitation data covers the period from
        $(input_meteoveg_reference_date + Day(startday_iso)) to $(input_meteoveg_reference_date + Day(stopday_iso))
        it will be cropped to the period determined by the other meteorologic inputs going from
        $(input_meteoveg_reference_date + Day(startday_veg)) to $(input_meteoveg_reference_date + Day(stopday_veg)).
        """
        # input_meteoiso = @chain input_meteoiso begin
        #    where(:days .>= startday, :days .<= endday) # Acutally not needed to crop
        #    end
        input_meteoveg_mod = input_meteoveg #Not needed to crop
    elseif ((stopday_iso < stopday_veg) | (startday_iso > startday_veg))
        @warn """
        Isotopic signature of precipitation data covers the period from
            $(input_meteoveg_reference_date + Day(startday_iso)) to $(input_meteoveg_reference_date + Day(stopday_iso))
        It does therefore not cover the entire period of other meteorologic inputs going from
            $(input_meteoveg_reference_date + Day(startday_veg)) to $(input_meteoveg_reference_date + Day(stopday_veg)).
        Simulation period will be limited to the period when isotopic signatures are available.
        Note that the initial conditions will be applied to the beginning of the new simulation period.
        """
        input_meteoveg_mod = filter(:days => (d) -> d .>= startday && d .<= endday, input_meteoveg)
    else
        input_meteoveg_mod = input_meteoveg #Not needed to crop
    end

    return input_meteoiso, input_meteoveg_mod
end

# path_initial_conditions = "examples/BEA2016-reset-FALSE-input/BEA2016-reset-FALSE_initial_conditions.csv"
function read_path_initial_conditions(path_initial_conditions)
    parsing_types =
        Dict(# Initial conditions (of vector states) -------
            "u_GWAT_init_mm" => Float64,       "u_INTS_init_mm" => Float64,
            "u_INTR_init_mm" => Float64,       "u_SNOW_init_mm" => Float64,
            "u_CC_init_MJ_per_m2"   => Float64,       "u_SNOWLQ_init_mm" => Float64)
    input_initial_conditions = DataFrame(File(path_initial_conditions;
        transpose=true, drop=[1], comment = "###",
        # Don't be strict, allow for NA as missing. Treat this later with disallowmissing!.
        types=parsing_types, missingstring = "NA"))

    expected_names = [String(k) for k in keys(parsing_types)]
    assert_colnames_as_expected(input_initial_conditions, path_initial_conditions, expected_names)

    # Impose type of Float64 instead of Float64?, by defining unused variables as -9999.99
    input_initial_conditions[2, "u_CC_init_MJ_per_m2"] = -9999.99
    input_initial_conditions[2, "u_SNOWLQ_init_mm"]    = -9999.99
    input_initial_conditions[3, "u_CC_init_MJ_per_m2"] = -9999.99
    input_initial_conditions[3, "u_SNOWLQ_init_mm"]    = -9999.99

    disallowmissing!(input_initial_conditions)

    return input_initial_conditions
end

# path_param = "examples/BEA2016-reset-FALSE-input/BEA2016-reset-FALSE_param.csv"
# TODO: check why read_path_param does not appear in the docs...
"""
    read_path_param(path_param; simulate_isotopes::Bool = false)

Reads in the `param.csv` based on a provided path. The `param.csv` has a structure shown in
the documentation (User Guide -> Input data). [Structure of input data](@ref)
"""
function read_path_param(path_param; simulate_isotopes::Bool = false)
    parsing_types =
        Dict(### Isotope transport parameters  -------,NA
            # "TODO" => Float64, "TODO2" => Float64,
            # TODO(bernhard): this needs to be extended with the currently hardcoded isotope transport parameters
            "DISPERSIVITY_mm" => Float64,
            "VXYLEM_mm" => Float64,
            # Meteorologic site parameters -------
            "LAT_DEG" => Float64,
            "ESLOPE_DEG" => Float64,       "ASPECT_DEG" => Float64,
            "ALB" => Float64,              "ALBSN" => Float64,
            "C1" => Float64,               "C2" => Float64,               "C3" => Float64,
            "WNDRAT" => Float64,           "FETCH" => Float64,            "Z0W" => Float64,              "ZW" => Float64,
            # Canopy parameters -------
            "MAXLAI" => Float64,
            "DENSEF_baseline_" => Float64,
            "SAI_baseline_" => Float64,
            "AGE_baseline_yrs" => Float64,
            "HEIGHT_baseline_m" => Float64,
            "LWIDTH" => Float64,           "Z0G" => Float64,              "Z0S" => Float64,
            "LPC" => Float64,              "CZS" => Float64,
            "CZR" => Float64,              "HS" => Float64,               "HR" => Float64,
            "ZMINH" => Float64,            "RHOTP" => Float64,            "NN" => Float64,
            # Interception parameters -------
            "FRINTLAI" => Float64,         "FSINTLAI" => Float64,
            "FRINTSAI" => Float64,         "FSINTSAI" => Float64,
            "CINTRL" => Float64,           "CINTRS" => Float64,
            "CINTSL" => Float64,           "CINTSS" => Float64,
            "RSTEMP" => Float64,
            # Snowpack parameters -------
            "MELFAC" => Float64,           "CCFAC" => Float64,            "LAIMLT" => Float64,
            "SAIMLT" => Float64,           "GRDMLT" => Float64,           "MAXLQF" => Float64,
            "KSNVP" => Float64,            "SNODEN" => Float64,
            # Leaf evaporation parameters (affecting PE) -------
            "GLMAX" => Float64,            "GLMIN" => Float64,            "CR" => Float64,               "RM" => Float64,
            "R5" => Float64,               "CVPD" => Float64,             "TL" => Float64,               "T1" => Float64,
            "T2" => Float64,               "TH" => Float64,
            # Plant parameters (affecting soil-water supply) -------
            "MXKPL" => Float64,
            "MXRTLN" => Float64,           "INITRLEN" => Float64,         "INITRDEP" => Float64,
            "RGRORATE" => Float64,         "RGROPER" => Float64,          "FXYLEM" => Float64,
            "PSICR" => Float64,            "RTRAD" => Float64,            "NOOUTF" => Int64,    # TODO(bernhard): make boolean
            # Soil parameters -------
            "IDEPTH_m" => Float64,           "QDEPTH_m" => Float64,
            "RSSA" => Float64,             "RSSB" => Float64,             "INFEXP" => Float64,
            "BYPAR" => Int64  ,             # TODO(bernhard): make boolean
            "QFPAR" => Float64,            "QFFC" => Float64,
            "IMPERV" => Float64,           "DSLOPE" => Float64,           "LENGTH_SLOPE" => Float64,
            "DRAIN" => Float64,            "GSC" => Float64,              "GSP" => Float64,
            # Numerical solver parameters -------
            "DTIMAX" => Float64,           "DSWMAX" => Float64,           "DPSIMAX" => Float64)
    # if (!simulate_isotopes)
    #     delete!(parsing_types, "TODO")
    #     delete!(parsing_types, "TODO2")
    # end

    input_param_df = DataFrame(File(path_param;
        transpose=true, drop=[1], comment = "###",
        # Be strict about loading NA's -> error if NA present
        types = parsing_types, missingstring = nothing, strict=true))

    expected_names = [String(k) for k in keys(parsing_types)]
    assert_colnames_as_expected(input_param_df, path_param, expected_names)

    # Set minimum/maximum values
    # from LWFBrook90R:PFILE.h
    input_param_df[:,:FXYLEM]  = min.(input_param_df[:,:FXYLEM], 0.990)
    input_param_df[:,:INITRLEN] = max.(input_param_df[:,:INITRLEN], 0.010)
    input_param_df[:,:INITRDEP] = max.(input_param_df[:,:INITRDEP], 0.010)

    # Impose type of Float64 instead of Float64?
    disallowmissing!(input_param_df)
    # Transform to NamedTuple
    input_param    = NamedTuple(input_param_df[1, Not([:DTIMAX, :DSWMAX, :DPSIMAX])])
    solver_opts    = NamedTuple(input_param_df[1, [:DTIMAX, :DSWMAX, :DPSIMAX]])

    return input_param, solver_opts
end

# path_storm_durations = "examples/BEA2016-reset-FALSE-input/BEA2016-reset-FALSE_meteo_storm_durations.csv"
function read_path_storm_durations(path_storm_durations)
    parsing_types =
        Dict("month" => String, "average_storm_duration_h" => Float64)
    input_storm_durations = DataFrame(File(path_storm_durations;
        comment = "###",
        # Be strict about loading NA's -> error if NA present
        types = parsing_types, missingstring = nothing, strict = true))

    expected_names = [String(k) for k in keys(parsing_types)]
    assert_colnames_as_expected(input_storm_durations, path_storm_durations, expected_names)

    # Assert months to be ordered
    received_month_names = input_storm_durations[!,"month"]
    expected_month_names = ["January", "Februrary", "March", "April", "May", "June", "July",
                            "August", "September", "October", "November", "December"]
    @assert received_month_names == expected_month_names """
        Input file '$path_storm_durations' does not contain the months in the expected order.
        Please correct the input file.
        \n\nExpected:\n$expected_month_names\nReceived:\n$received_month_names\n
        """

    # Rename column names
    rename!(input_storm_durations,
        :average_storm_duration_h => :storm_durations_h)

    return input_storm_durations
end
# path_soil_horizons = "examples/BEA2016-reset-FALSE-input/BEA2016-reset-FALSE_soil_horizons.csv"
function read_path_soil_horizons(path_soil_horizons)
    # Derive whether to use Clapp Hornberger or MualemVanGenuchten based on the input data
    MualVanGen_expected_column_names =
        ["HorizonNr","Upper_m","Lower_m","ths_volFrac","thr_volFrac","alpha_perMeter","npar_","ksat_mmDay","tort_","gravel_volFrac"]
    ClappHornberger_expected_column_names =
        ["HorizonNr","Upper_m","Lower_m","thsat_volFrac","thetaf_volFrac","psif_kPa","bexp_","kf_mmDay","wtinf_","gravel_volFrac"]

    if strip.(String.(propertynames(File(path_soil_horizons)))) == MualVanGen_expected_column_names
        FLAG_MualVanGen = 1
    elseif strip.(String.(propertynames(File(path_soil_horizons)))) == ClappHornberger_expected_column_names
        FLAG_MualVanGen = 0
    else
        error("""
        Could not derive which hydraulic model parametrization (Mualem-Van-Genuchten or
        Clapp-Hornberger) to use, based the column names of the input file '$path_soil_horizons'.
        Please check and correct the input file!

        Expected column names are either:
        Mualem-Van-Genuchten: $MualVanGen_expected_column_names
        or for Clapp-Hornberger: $ClappHornberger_expected_column_names
        """)
    end

    if FLAG_MualVanGen == 1
        # Assert units:
        assert_unitsHeader_as_expected(path_soil_horizons,
            DataFrame(HorizonNr = "-", Upper_m = "m", Lower_m = "m", ths_volFrac = "volume fraction (-)",
                thr_volFrac = "volume fraction (-)", alpha_perMeter = "perMeter", npar_ = "-",
                ksat_mmDay = "mm per day", tort_ = "-", gravel_volFrac = "volume fraction (-)"))

        # Load file
        input_soil_horizons = DataFrame(File(path_soil_horizons;
            skipto=3, header=MualVanGen_expected_column_names,
            delim=',',
            # Be strict about loading NA's -> error if NA present
            types=[Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64],
            missingstring = nothing, strict = true))

    else # FLAG_MualVanGen = 0 (Clapp-Hornberger)
        # Assert units:
        assert_unitsHeader_as_expected(path_soil_horizons,
            DataFram(HorizonNr = "-", Upper_m = "m", Lower_m = "m", thsat_volFrac = "volume fraction (-)",
                thetaf_volFrac = "volume fraction (-)", psif_kPa = "kPa", bexp_ = "-",
                kf_mmDay = "mm per day", wtinf_ = "-", gravel_volFrac = "volume fraction "))

        # Load file
        input_soil_horizons = DataFrame(File(path_soil_horizons;
            skipto=3, header=ClappHornberger_expected_column_names,
            delim=',', # ignorerepeated=true
            # Be strict about loading NA's -> error if NA present
            types=[Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64],
            missingstring = nothing, strict = true))
    end

    # Check that defined horizons do not overlap
    @assert input_soil_horizons[1:end-1,"Lower_m"] == input_soil_horizons[2:end,"Upper_m"] """
        Input file '$path_soil_horizons' contains overlapping layers.
        Please check and correct the input file.
        """

    @assert input_soil_horizons[1,"Upper_m"] ≈ 0 """
        Input file '$path_soil_horizons' does not start with 0.0 as Upper_m limit of first horizon.
        Please check and correct the input file.
        """

    # Assert that no missing values in
    # Impose type of Float64 instead of Float64?
    disallowmissing!(input_soil_horizons)

    # Make dataframe representing physical horizons/layers in 1D domain
    # soil_horizons = DataFrame()
    soil_horizons = input_soil_horizons[:,["HorizonNr", "Upper_m", "Lower_m"]];
    soil_horizons[!, :shp] = LWFBrook90.MualemVanGenuchtenSHP(input_soil_horizons);

    return soil_horizons
end

# path_soil_discretization = "examples/BEA2016-reset-FALSE-input/BEA2016-reset-FALSE_soil_discretization.csv"
function read_path_soil_discretization(path_soil_discretization)
    parsing_types =
        Dict("Upper_m"      => Float64,
             "Lower_m"      => Float64,
             "Rootden_"      => Float64,
             "uAux_PSIM_init_kPa"   => Float64,
             "u_delta18O_init_permil" => Float64,
             "u_delta2H_init_permil"  => Float64)

    input_soil_discretization = DataFrame(File(path_soil_discretization;
        skipto=3, missingstring = "NA", types=parsing_types))

    # Assert colnames
    expected_names = [String(k) for k in keys(parsing_types)]
    assert_colnames_as_expected(input_soil_discretization, path_soil_discretization, expected_names)

    # Check that defined layers do not overlap
    @assert input_soil_discretization[1:end-1,"Lower_m"] == input_soil_discretization[2:end,"Upper_m"] """
        Input file '$path_soil_discretization' contains overlapping layers.
        Please check and correct the input file.
        """
    @assert input_soil_discretization[1,"Upper_m"] ≈ 0 """
        Input file '$path_soil_discretization' does not start with 0.0 as Upper_m limit of first layer.
        Please check and correct the input file.
        """
    # Check that defined layers are not zero height
    zero_height_layers = input_soil_discretization[1:end-1,"Upper_m"] .== input_soil_discretization[2:end,"Upper_m"]
    idx_layer_zero_height = (1:nrow(input_soil_discretization)-1)[zero_height_layers]
    @assert !any(zero_height_layers) """
        Input file '$path_soil_discretization' contains layers of zero height:
        $(input_soil_discretization[[idx_layer_zero_height[1]-1, idx_layer_zero_height[1], idx_layer_zero_height[1]+1], :])
        Please check and correct the input file.
        """

    # Assert units:
    assert_unitsHeader_as_expected(path_soil_discretization,
        DataFrame(Upper_m = "m", Lower_m = "m", Rootden_ = "-",
            uAux_PSIM_init_kPa = "kPa",
            u_delta18O_init_permil = "permil", u_delta2H_init_permil = "permil"))

    return input_soil_discretization
end

function assert_colnames_as_expected(input_df, input_path, expected_names)
    # Assert that inputs are correct:
    received_names = names(input_df)

    @assert all(expected_names .∈ (received_names,)) """
    Input file '$input_path' does not contain all the expected column names.
    \nExpected but missing:\n$(expected_names[(!).(expected_names .∈ (received_names,))])
    \n\nExpected:\n$expected_names\nReceived:\n$received_names\n
    """
    @assert all(received_names .∈ (expected_names,)) """
    Input file '$input_path' contains other than the expected  column names.
    \nReceived unexpected: $(received_names[(!).(received_names .∈ (expected_names,))])
    \nExpected:\n$expected_names\nReceived:\n$received_names\n
    """
end

function assert_unitsHeader_as_expected(path, expected_units)
    actual_units = DataFrame(File(path;skipto=2,limit=1,types=String))
    @assert actual_units == expected_units """
    Unexpected units in input file $path. Expected units:
    $expected_units
    Received units:
    $actual_units"""
end

function assert_unitsColumn_as_expected(path, expected_units)
    # TODO(bernhard): define units in and implement this check.
    # @assert DataFrame(File(path;skipto=2,limit=1)) == expected_units """
    # Unexpected units in input file $path. Expected units:
    # $expected_units"""
end
