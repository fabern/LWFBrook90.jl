using CSV: read, File
using DataFrames: DataFrame, rename, sort!# ,select
using DataFramesMeta#: @linq, transform, DataFramesMeta
using Dates: DateTime, Millisecond, Second, Day, Month, month, value, dayofyear, format
using Statistics: mean
# using Printf: @sprintf
"""
    loadSPAC(folder::String, prefix::String)

Define instance of SPAC model by loading different input files for LWFBrook90:
- `meteoveg.csv`
- `param.csv`
- `meteo_storm_durations.csv`
- `initial_conditions.csv`
- `soil_horizons.csv`

These files were created with an R script `generate_LWFBrook90jl_Input.R` that
takes the same arguements as the R function `LWFBrook90R::run_LWFB90()` and generates
the corresponding input files for LWFBrook90.jl.
"""
function loadSPAC(folder::String, prefix::String;
    suffix::String = "",
    simulate_isotopes::Bool = true,
    compute_intermediate_quantities::Bool = true)

    ## Define solver options
    solver_options =
        (#Reset                           = false, # currently only Reset = 0 implemented
         compute_intermediate_quantities = true,   # Flag whether ODE containes additional quantities than only states
         simulate_isotopes               = simulate_isotopes
        )

    ## Define paths of all input files
    input_file_XXXX = prefix*"_XXXX"*suffix*".csv"
    path_meteoveg           = joinpath(folder, replace(input_file_XXXX, "XXXX" => "meteoveg"))
    path_param              = joinpath(folder, replace(input_file_XXXX, "XXXX" => "param"))
    # unused path_precdat   = joinpath(folder, replace(input_file_XXXX, "XXXX" => "precdat"))
    path_storm_durations    = joinpath(folder, replace(input_file_XXXX, "XXXX" => "meteo_storm_durations"))
    path_initial_conditions = joinpath(folder, replace(input_file_XXXX, "XXXX" => "initial_conditions"))
    path_soil_horizons      = joinpath(folder, replace(input_file_XXXX, "XXXX" => "soil_horizons"))

    ## Load time-varying atmospheric forcing
    reference_date, tspan, input_meteoveg, meteo_iso_forcing, storm_durations =
        init_forcing(path_meteoveg, path_storm_durations; simulate_isotopes)

    meteo_forcing = input_meteoveg[:, [:days, :GLOBRAD, :TMAX, :TMIN, :VAPPRES, :WIND, :PRECIN]]

    ## Load time-varying vegetation parameters
    if ("DENSEF" in names(input_meteoveg) &&
        "HEIGHT" in names(input_meteoveg) &&
        "LAI" in names(input_meteoveg) &&
        "SAI" in names(input_meteoveg))
        _canopy_evolution = input_meteoveg[:, [:days, :DENSEF, :HEIGHT, :LAI, :SAI]]
    else
        # This case is currently not implemented
        @warn  """
                Input_meteoveg is expected to contain one or multiple of the columns: :DENSEF, :HEIGHT, :LAI, or :SAI.
                Please check your input files with the current documentation and possibly contact the developer team if the error persists.
                If it does not contain the columns please provide a parametrization to loadSPAC(, canopy_evolution::NamedTuple = ())
                with the named tuple:
                `(DENSEF = 100, HEIGHT = 25, SAI = 100, LAI = (DOY_Bstart, Bduration, DOY_Cstart, Cduration, LAI_perc_BtoC, LAI_perc_CtoB))`
                """
        @assert canopy_evolution isa NamedTuple
        @assert keys(canopy_evolution) == (:DENSEF, :HEIGHT, :SAI, :LAI)
        @assert keys(canopy_evolution.LAI) == (:DOY_Bstart, :Bduration, :DOY_Cstart, :Cduration, :LAI_perc_BtoC, :LAI_perc_CtoB)
        _canopy_evolution = canopy_evolution
        # canopy_evolution = ... # Will be overwritten with a dataframe, when discretizing loadSPAC()
        #       TODO in setup(): DataFrame(days=input_meteoveg.days, DENSEF=100, HEIGHT=100, LAI = NaN, SAI=100))
        #       TODO in setup(): canopy_evolution.LAI = interpolate_LAI(input_meteoveg.days, LAI_DOY_B1, LAI_DOY_B2, LAI_DOY_C1, LAI_DOY_C2, LAImin_percent, LAImax_percent) #TODO: make this vector by interpolating between input_meteoveg.days as DOY

        ### expected format of `canopy_evolution.csv`:
        ### densef_percent, height_percent, SAI_percent, LAImin_percent, LAI_DOY1_budburst, LAI_duration1, LAI_DOY2_leaffall, LAI_duration2)
        ### 100,            100,            100,         20,             130,               160,          270,               320

        # LAI_perc_BtoC
        # LAI_perc_CtoB
        # LAI_DOY_B1 = DOY_Bstart
        # LAI_DOY_B2 = DOY_Bstart + Bduration
        # LAI_DOY_C1 = DOY_Cstart
        # LAI_DOY_C2 = DOY_Cstart + Cduration
        # @assert 1 ≤ LAI_DOY_B1 ≤ LAI_DOY_B2 ≤ LAI_DOY_C1 ≤ LAI_DOY_C2 ≤ 365
            # example R-code:
            # compute_LAI_percent <- function(DOY, LAI_min, DOYinit, DOY1, DOY2, DOY3, DOY4, DOYmax){
            #   LAI_max <- 1
            #   LAI_dates  <- c(DOYinit,    DOY1, DOY2   , DOY3   , DOY4   , DOYmax)
            #   LAI_values <- c(LAI_min, LAI_min, LAI_max, LAI_max, LAI_min, LAI_min)
            #   # return relative evolution of LAI
            #   100*stats::approx(x = LAI_dates, y = LAI_values, method = "linear", xout = DOY)$y
            # }
    end


    ## Load space-varying soil data
    soil_horizons = init_soil(path_soil_horizons)

    ## Define placeholder for root distribution (to be loaded either via file or then parametrized)
    root_distribution = (input = "input_soildiscretization",
                         distribution = DataFrame())
    root_distribution = (input = (rootden = ((Δz_m)->LWFBrook90.Rootden_beta_(0.97, Δz_m = Δz_m)),
                         distribution = DataFrame()))
                            # Δz = (thickness_m = nothing,
                            #         functions   = (rootden   = ((Δz_m)->LWFBrook90.Rootden_beta_(0.97, Δz_m = Δz_m)),
                            #                 PSIM_init = ((Δz_m)->fill(-6.3, length(Δz_m))),
                            #                 δ18Ο_init = ((Δz_m)->ifelse.(cumsum(Δz_m) .<= 0.2, -13., -10.)),
                            #                 δ2Η_init  = ((Δz_m)->ifelse.(cumsum(Δz_m) .<= 0.2, -95., -70.))))

    ## Load initial conditions of scalar state variables
    IC = init_IC(path_initial_conditions; simulate_isotopes)

    ## Load model input parameters
    params, solver_opts = init_param(path_param; simulate_isotopes = simulate_isotopes)

    solver_options = merge(solver_options, solver_opts) # append to manually provided solver options

    ## Make time dependent input parameters continuous in time (interpolate)
    (meteo_forcing_cont, meteo_iso_forcing_cont) =
        LWFBrook90.interpolate_meteo(;
            meteo_forcing                 = meteo_forcing,
            meteo_iso_forcing             = meteo_iso_forcing);

    return SPAC(;
        reference_date    = reference_date,
        tspan             = tspan,
        solver_options    = solver_options,
        Δz_thickness_m    = IC.soil, # path to soil_discretization.csv
                            # - `Δz_thickness_m`:either
                            #     - `String` defining the path to a soil_discretizations.csv
                            #     - `DataFrame` of contents of soil_discretizations.csv
                            #     - or `Vector`: Δz_thickness_m
        forcing = (meteo           = meteo_forcing_cont,
                   meteo_iso       = meteo_iso_forcing_cont,
                   storm_durations = storm_durations),
        pars    = (params = params,
                   root_distribution = root_distribution,
                   IC_scalar = IC.scalar,
                   IC_soil = "Δz_thickness_m",
                   canopy_evolution = _canopy_evolution,
                   soil_horizons = soil_horizons),
        )
end
# -   `pars`: NamedTuple: (:params, )
#     - `pars.params`: NamedTuple: (), containing scalar parameter values for the model
#     - `pars.root_distribution`: either:
#         - NamedTuple: (beta = 0.97, z_rootMax_m = nothing) parametrization of root distribution with depth (f(z_m) with z_m in meters and negative downward). Alternatively path to soil_discretization.csv" #TODO: interpolate this in time
#         - or String: "Δz_thickness_m" meaning that it must be defined in soil_discretizations.csv
#     - `pars.IC_scalar`: NamedTuple: (), containing initial conditions of the scalar state variables
#     - `pars.IC_soil`: initial conditions of the state variables (scalar or related to the soil). Either:
#         - NamedTuple: (PSIM_init, δ18O_init, δ2H_init), containing the initial values
#         - NamedTuple: (PSIM_init_fct, δ18O_init_fct, δ2H_init_fct), containing the initial values
#         - or String: "Δz_thickness_m" meaning that it must be defined in soil_discretizations.csv
#     - `pars.canopy`: canopy parameters (LAI, SAI, DENSEF, HEIGHT). Either:
#         - NamedTuple: (DENSEF = 100, HEIGHT = 25, SAI = 100, LAI = (DOY_Bstart, Bduration, DOY_Cstart, Cduration, LAI_perc_BtoC, LAI_perc_CtoB)), containing constant values and parameters for LAI_relative interpolation
#         - or String: "Δz_thickness_m" meaning that it must be defined in soil_discretizations.csv
#     - `pars.soil_horizons`: DataFrame containing description of soil layers/horizons and soil hydraulic parameters
# """
# SPAC(reference_date,
#     tspan,
#     solver_options,                 # reset_flag, compute_intermediate_quantities, simulate_isotopes, soil_output_depths
#     Δz_thickness_m, # either "path_to_soil_discretization.csv" or then NamedTuple: `(Δz_thickness_m = , functions = )`
#     # A2) Physical forcing
#     forcing, # atmospheric forcing:  #TODO: interpolate forcing.meteo and forcing.meteo_iso in time
#     # B) Model parameters: might be estimated
#     pars,
#         # pars.params: e.g. (VXYLEM_mm = 20.0, DISPERSIVITY_mm = 10, Str = "oh")
#         # pars.root_distribution: e.g. (beta = 0.97, z_rootMax_m = nothing)
#         # pars.IC:                e.g. .soil = PSIM, δ2H, δ18O, but also .scalar = SNOW,INTR,INTS
#         # pars.root_distribution = (beta = 0.97, z_rootMax_m = nothing) # TODO: define default....
#         # pars.canopy:            e.g. LAI(t), ...
#         # pars.soil_horizons
# )


function Base.show(io::IO, mime::MIME"text/plain", model::SPAC)
    println(io, "SPAC model:")
    println(io, "===== DATES:===============")
    println(io, format(model.reference_date, "YYYY-mm-dd"), ", tspan:", model.tspan, "days")
    println(io, "\n===== METEO FORCING:===============")
    show_avg_and_range = function(vector, title)
        # "$(title)avg:$(round(Statistics.mean(vector); digits=2)), range:$(extrema(vector))"
        @sprintf("%savg:%7.2f, range:%7.2f to%7.2f",
                 title, mean(vector), extrema(vector)[1], extrema(vector)[2])
    end
    println(io, show_avg_and_range(model.forcing.meteo.p_GLOBRAD.itp.itp.coefs,     "GLOBRAD (MJ/m2/day): "))
    println(io, show_avg_and_range(model.forcing.meteo.p_PREC.itp.itp.coefs,        "PREC       (mm/day): "))
    println(io, show_avg_and_range(model.forcing.meteo.p_TMAX.itp.itp.coefs,        "TMAX           (°C): "))
    println(io, show_avg_and_range(model.forcing.meteo.p_TMIN.itp.itp.coefs,        "TMIN           (°C): "))
    println(io, show_avg_and_range(model.forcing.meteo.p_VAPPRES.itp.itp.coefs,     "VAPPRES       (kPa): "))
    println(io, show_avg_and_range(model.forcing.meteo.p_WIND.itp.itp.coefs,        "WIND          (m/s): "))
    if (model.solver_options.simulate_isotopes)
        println(io, show_avg_and_range(model.forcing.meteo_iso.p_d18OPREC.itp.coefs,"δ18O            (‰): "))
        println(io, show_avg_and_range(model.forcing.meteo_iso.p_d2HPREC.itp.coefs, "δ2H             (‰): "))
    end
    # println(io, model.forcing.storm_durations)

    println(io, "\n===== CANOPY EVOLUTION:===============")
    if model.pars.canopy_evolution isa DataFrame
        println(io, "model.pars.canopy_evolution was loaded form meteoveg.csv")
        # println(io, show_avg_and_range(model.pars.canopy_evolution.DENSEF,"DENSEF            (%): "))
        # println(io, show_avg_and_range(model.pars.canopy_evolution.HEIGHT,"HEIGHT            (%): "))
        # println(io, show_avg_and_range(model.pars.canopy_evolution.LAI,"LAI            (%): "))
        # println(io, show_avg_and_range(model.pars.canopy_evolution.SAI,"SAI            (%): "))
    else
        println(io, model.pars.canopy_evolution)
    end

    println(io, "\n===== INITIAL CONDITIONS:===============")
    println(io, "Soil   IC: " * model.pars.IC_soil)
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
    print(  io, "Δz_thickness_m:    "); println(io, model.Δz_thickness_m)
    print(  io, "Root distribution: "); println(io, model.pars.root_distribution)
    println(io, "Soil layer properties:")
    # println(io, model.pars.soil_horizons)
    # for shp in model.pars.soil_horizons.shp println(io, shp) end
    # show(IOContext(io, :limit => false), "text/plain", model.pars.soil_horizons)
    show(IOContext(io, :limit => false), mime, model.pars.soil_horizons[:,Not(:HorizonNr)]; truncate = 100);
    println(io, "")

    println(io, "\n===== SOLVER OPTIONS:===============")
    println(io, model.solver_options)
end

######################
# Define functions to handle DateTimes and convert into Days as Floats
"""
    DateTime2RelativeDaysFloat(x,reference_DateTime)

Transforms DateTimes `x` to simulation time
"""
function DateTime2RelativeDaysFloat(x::DateTime, reference::DateTime)
    ms2days = 1.0/(24.0*3600.0*1000.0) # to convert milliseconds to days
    ms2days*value(convert(Millisecond, x-reference))
end
"""
    RelativeDaysFloat2DateTime(t, reference_DateTime)

Transforms simulation time `t` to DateTimes
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

function init_forcing(path_meteoveg, path_storm_durations; simulate_isotopes = true)
    if (simulate_isotopes)
        path_meteoiso = replace(path_meteoveg, "meteoveg" => "meteoiso")
    end
    # if (simulate_isotopes && prefix == "2021-09-23-DAV")
    #     # TODO(bernhard): hardcoded overriding of Davos example data
    #     path_meteoiso = joinpath("2021-09-23_onlyDavos2020/2021-09-23-DAV-isoInput", "2021-10-04-DAV_meteoIso"*suffix*".csv")
    #         path_initial_conditions = joinpath("2021-09-23_onlyDavos2020/2021-09-23-DAV-isoInput",
    #                                             "2021-09-23-DAV_initial_conditions_tracers.csv")
    #         path_param              = joinpath("2021-09-23_onlyDavos2020/2021-09-23-DAV-isoInput",
    #                                             "2021-09-23-DAV_paramIso.csv")
    # end


    meteo_forcing, reference_date = read_path_meteoveg(path_meteoveg)

    if (simulate_isotopes)
        meteo_iso_forcing, meteo_forcing = read_path_meteoiso(path_meteoiso,
            meteo_forcing,
            reference_date)
    else
        meteo_iso_forcing = nothing
    end

    tspan = extrema(meteo_forcing.days)
    # tspan = convert.(Int64, extrema(meteo_forcing.days))

    # Load precipitation subdaily details
    #' @param storm_durations a [1,12]-matrix of precipitation durations (hours) for each month.
    storm_durations = read_path_storm_durations(path_storm_durations)
    #' @param precdat A matrix of precipitation interval data with 5 columns:
    #'                   year, month, day, interval-number (1:precint), prec.
    # input_precdat = DataFrame(a = Nothing, b = Nothing)
    # TODO(bernhard): currently not implemented.
    #                 Only using PRECIN (resolution daily).
    #                 PRECDAT would allow to have smaller resolution (would require changes).
    # unused: input_precdat = read(path_precdat, DataFrame;
    # unused:           missingstring = "-999",
    # unused:           skipto=2, header=["year","month","day","interval_number","prec"], delim=',',
    # unused:           ignorerepeated=true)

    # Assert that no missing
    # Impose type of Float64 instead of Float64?
    disallowmissing!(meteo_forcing)
    if (simulate_isotopes)
        disallowmissing!(meteo_iso_forcing, [:days, :delta18O_permil, :delta2H_permil])
    end
    disallowmissing!(storm_durations, [:month, :storm_durations_h])

    return reference_date, tspan, meteo_forcing, meteo_iso_forcing, storm_durations
end

function init_IC(path_initial_conditions; simulate_isotopes = true)
    scalar_initial_conditions = read_path_initial_conditions(path_initial_conditions;
                                                                simulate_isotopes = simulate_isotopes)
    # Assert that no missing
    # Impose type of Float64 instead of Float64?
    disallowmissing!(scalar_initial_conditions)

    # Note: IC.soil is used this if there is a definition of the initial conditions that is independent from the discretization
    #       If there is no such definition of the ICs, set soil to the path of the soil_discretiztation file that will be
    #       read definition of soil initial condition (uAux_PSIM_init_kPa, as well as u_delta18O_init_permil, u_delta2H_init_permil)
    # path_soil_discretization = joinpath(folder, replace(input_file_XXXX, "XXXX" => "soil_discretization"))
    path_soil_discretization = replace(path_initial_conditions, "initial_conditions" => "soil_discretization")

    IC = (soil   = path_soil_discretization,    # if set to a path, this will be loaded when discretizing the soil domain
          scalar = scalar_initial_conditions)
    return IC
end

function init_soil(path_soil_horizons)
    read_path_soil_horizons(path_soil_horizons)
end

function init_param(path_param; simulate_isotopes = true)
    input_param, solver_opts = read_path_param(path_param; simulate_isotopes = simulate_isotopes)
end

# path_meteoveg = "examples/BEA2016-reset-FALSE-input/BEA2016-reset-FALSE_meteoveg.csv"
function read_path_meteoveg(path_meteoveg)
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

    input_meteoveg = @linq DataFrame(File(path_meteoveg;
        skipto=3, delim=',', ignorerepeated=false,
        # Be strict about loading NA's -> error if NA present
        types=parsing_types, missingstring = nothing, strict=true))  |>
        transform(:dates = DateTime.(:dates))

    expected_names = [String(k) for k in keys(parsing_types)]
    assert_colnames_as_expected(input_meteoveg, path_meteoveg, expected_names)

    # Assert units:
    assert_unitsHeader_as_expected(path_meteoveg,
        DataFrame(dates = "YYYY-MM-DD", globrad_MJDayM2 = "MJ/Day/m2",
        tmax_degC = "degree C", tmin_degC = "degree C", vappres_kPa = "kPa",
        windspeed_ms = "m per s", prec_mmDay = "mm per day",
        densef_percent = "percent", height_percent = "percent",
        lai_percent = "percent", sai_percent = "percent"))

    # Assert validity of values
    @assert all(input_meteoveg.densef_percent .> 5) """
        Densef_percent in meteoveg.csv should not be set lower than 5% as it affects aerodynamics.
    """

    # Assert that no gaps
    @assert all(diff(input_meteoveg.dates) .== Millisecond(86400000)) "There are gaps in the forcing file. The file ($path_meteoveg) needs to have a value for each data from start until the end."

    # Identify period of interest
    # Starting date: latest among the input data
    # Stopping date: earliest among the input data
    starting_date = maximum(minimum,[input_meteoveg[:,"dates"]])
    stopping_date = minimum(maximum,[input_meteoveg[:,"dates"]])

    input_meteoveg = @linq input_meteoveg |>
        where(:dates .>= starting_date, :dates .<= stopping_date)

    rename!(input_meteoveg,
        :globrad_MJDayM2 => :GLOBRAD,
        :tmax_degC       => :TMAX,
        :tmin_degC       => :TMIN,
        :vappres_kPa     => :VAPPRES,
        :windspeed_ms    => :WIND,
        :prec_mmDay      => :PRECIN,
        :densef_percent  => :DENSEF,
        :height_percent  => :HEIGHT,
        :lai_percent     => :LAI,
        :sai_percent     => :SAI)

    # Transform times from DateTimes to simulation time (Float of Days)
    reference_date = starting_date

    input_meteoveg = @linq input_meteoveg |>
        transform(:dates = DateTime2RelativeDaysFloat.(:dates, reference_date)) |>
        rename(Dict(:dates => :days))

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

        input_meteoiso = @linq DataFrame(File(path_meteoiso; header = 4,
                    skipto=5, delim=',', ignorerepeated=true,
                    # Don't be strict, allow for NA as missing. Treat this later with disallowmissing!.
                    types=parsing_types, missingstring = ["","NA"]))
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
        input_meteoiso = @linq input_meteoiso |>
                transform(:dates = ceil.(:dates, Month))
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
    input_meteoiso = @linq input_meteoiso |>
        transform(:dates = DateTime2RelativeDaysFloat.(:dates, input_meteoveg_reference_date)) |>
        rename(Dict(:dates => :days))

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
        # input_meteoiso = @linq input_meteoiso |>
        #    where(:days .>= startday, :days .<= endday) # Acutally not needed to crop
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
        input_meteoveg_mod = @linq input_meteoveg |>
            where(:days .>= startday, :days .<= endday)
    else
        input_meteoveg_mod = input_meteoveg #Not needed to crop
    end

    return input_meteoiso, input_meteoveg_mod
end

# path_initial_conditions = "examples/BEA2016-reset-FALSE-input/BEA2016-reset-FALSE_initial_conditions.csv"
function read_path_initial_conditions(path_initial_conditions; simulate_isotopes::Bool = false)
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


    if String.(propertynames(File(path_soil_horizons))) == MualVanGen_expected_column_names
        FLAG_MualVanGen = 1
    elseif String.(propertynames(File(path_soil_horizons))) == ClappHornberger_expected_column_names
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
    @assert DataFrame(File(path;skipto=2,limit=1)) == expected_units """
    Unexpected units in input file $path. Expected units:
    $expected_units"""
end

function assert_unitsColumn_as_expected(path, expected_units)
    # TODO(bernhard): define units in and implement this check.
    # @assert DataFrame(File(path;skipto=2,limit=1)) == expected_units """
    # Unexpected units in input file $path. Expected units:
    # $expected_units"""
end
