using CSV: read, File
using DataFrames: DataFrame, rename, sort!# ,select
using DataFramesMeta#: @linq, transform, DataFramesMeta
using Dates: DateTime, Millisecond, Second, Day, Month, month, value, dayofyear

"""
    read_inputData(folder::String, prefix::String)

Load different input files for LWFBrook90:
- meteoveg.csv
- param.csv
- meteo_storm_durations.csv
- initial_conditions.csv
- soil_horizons.csv
- soil_discretization.csv

These files were created with an R script `generate_LWFBrook90Julia_Input.R` that
takes the same arguements as the R function `LWFBrook90R::run_LWFB90()` and generates
the corresponding input files for LWFBrook90.jl.
"""
# folder = input_path
# prefix = input_prefix
# suffix = ""
function read_inputData(folder::String, prefix::String;
    suffix::String = "",
    simulate_isotopes::Bool = false)
    # https://towardsdatascience.com/read-csv-to-data-frame-in-julia-programming-lang-77f3d0081c14
    ## A) Define paths of all input files
    path_meteoveg           = joinpath(folder, prefix*"_meteoveg"*suffix*".csv")
    path_param              = joinpath(folder, prefix*"_param"*suffix*".csv")
    # unused path_precdat=joinpath(folder, prefix*"_precdat"*suffix*".csv")
    path_storm_durations    = joinpath(folder, prefix*"_meteo_storm_durations"*suffix*".csv")
    path_initial_conditions = joinpath(folder, prefix*"_initial_conditions"*suffix*".csv")
    path_soil_horizons      = joinpath(folder, prefix*"_soil_horizons"*suffix*".csv")

    if (simulate_isotopes)
        path_meteoiso = joinpath(folder, prefix*"_meteoiso"*suffix*".csv")
    end
    if (simulate_isotopes && prefix == "2021-09-23-DAV")
        # TODO(bernhard): hardcoded overriding of Davos example data
        path_meteoiso = joinpath("2021-09-23_onlyDavos2020/2021-09-23-DAV-isoInput", "2021-10-04-DAV_meteoIso"*suffix*".csv")
            path_initial_conditions = joinpath("2021-09-23_onlyDavos2020/2021-09-23-DAV-isoInput",
                                                "2021-09-23-DAV_initial_conditions_tracers.csv")
            path_param              = joinpath("2021-09-23_onlyDavos2020/2021-09-23-DAV-isoInput",
                                                "2021-09-23-DAV_paramIso.csv")
    end

    ## B) Load input data (time- and/or space-varying parameters)
    input_meteoveg, input_meteoveg_reference_date = read_path_meteoveg(path_meteoveg)

    if (simulate_isotopes)
        input_meteoiso, input_meteoveg = read_path_meteoiso(path_meteoiso,
            input_meteoveg,
            input_meteoveg_reference_date)
    else
        input_meteoiso = nothing
    end

    ## C) Load other parameters
    # Load model input parameters
    #' @param param A numeric vector of model input parameters. Order:
    input_param = read_path_param(path_param; simulate_isotopes = simulate_isotopes)

    # Load initial conditions of scalar state variables
    input_initial_conditions = read_path_initial_conditions(path_initial_conditions)

    # Load soil data
    input_soil_horizons, simOption_FLAG_MualVanGen = read_path_soil_horizons(path_soil_horizons)

    # Load precipitation subdaily details
    #' @param storm_durations a [1,12]-matrix of precipitation durations (hours) for each month.
    input_storm_durations = read_path_storm_durations(path_storm_durations)
    #' @param precdat A matrix of precipitation interval data with 5 columns:
    #'                   year, month, day, interval-number (1:precint), prec.
    # input_precdat = DataFrame(a = Nothing, b = Nothing)
    # TODO(bernhard): currently not implemented.
    #                 Only using PRECIN (resolution daily).
    #                 PRECDAT would allow to have smaller resolution (would require changes).
    # unused: input_precdat = read(path_precdat, DataFrame;
    # unused:           missingstring = "-999",
    # unused:           datarow=2, header=["year","month","day","interval_number","prec"], delim=',',
    # unused:           ignorerepeated=true)

    return (input_meteoveg,
        input_meteoiso,
        input_meteoveg_reference_date,
        input_param,
        input_storm_durations,
        input_initial_conditions,
        input_soil_horizons,
        simOption_FLAG_MualVanGen)
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

# Subset input data and transform dates into floats relative to reference_date
# """
#    subset_standardize(data::DataFrame, start::DateTime, stop::DateTime, reference::DateTime)
#
# Returns DataFrame `data` that is subset between `start` and `stop` and colum `dates` transformed to simulation time.
# """
# function subset_standardize(data::DataFrame, start::DateTime, stop::DateTime, reference::DateTime)
#     @linq data |>
#     # Subset
#     where(:dates .>= start, :dates .<= stop) |>
#     # Compute time relative to reference date
#     transform(dates = DateTime2RelativeDaysFloat.(:dates, reference)) |>
#     # Rename colum
#     rename(Dict(:dates => :days))
# end

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
            :densef_         => Float64,
            :height_m        => Float64,
            :lai_            => Float64,
            :sai_            => Float64,
            :age_yrs         => Float64)

    input_meteoveg = @linq DataFrame(File(path_meteoveg;
        datarow=3, delim=',', ignorerepeated=true, types=parsing_types))  |>
        transform(dates = DateTime.(:dates))

    expected_names = [String(k) for k in keys(parsing_types)]
    assert_colnames_as_expected(input_meteoveg, path_meteoveg, expected_names)

    rename!(input_meteoveg,
        :globrad_MJDayM2 => :GLOBRAD,
        :tmax_degC       => :TMAX,
        :tmin_degC       => :TMIN,
        :vappres_kPa     => :VAPPRES,
        :windspeed_ms    => :WIND,
        :prec_mmDay      => :PRECIN,
        :densef_         => :DENSEF,
        :height_m        => :HEIGHT,
        :lai_            => :LAI,
        :sai_            => :SAI,
        :age_yrs         => :AGE)

    # Assert units:
    assert_unitsHeader_as_expected(path_meteoveg,
        DataFrame(dates = "YYYY-MM-DD", globrad_MJDayM2 = "MJ/Day/m2",
        tmax_degC = "degree C", tmin_degC = "degree C", vappres_kPa = "kPa",
        windspeed_ms = "m per s", prec_mmDay = "mm per day", densef_ = "-",
        height_m = "m", lai_ = "-", sai_ = "-", age_yrs = "years"))

    # Identify period of interest
    # Starting date: latest among the input data
    # Stopping date: earliest among the input data
    starting_date = maximum(minimum,[input_meteoveg[:,"dates"]])
    stopping_date = minimum(maximum,[input_meteoveg[:,"dates"]])

    input_meteoveg = @linq input_meteoveg |>
        where(:dates .>= starting_date, :dates .<= stopping_date)

    # Transform times from DateTimes to simulation time (Float of Days)
    reference_date = starting_date

    input_meteoveg = @linq input_meteoveg |>
        transform(dates = DateTime2RelativeDaysFloat.(:dates, reference_date)) |>
        rename(Dict(:dates => :days))

    return input_meteoveg, reference_date
end

function read_path_meteoiso(path_meteoiso,
        input_meteoveg,
        input_meteoveg_reference_date)

    parsing_types =
            Dict(:Site           => Char,
                :Date            => DateTime,
                :Latitude        => Float64,
                :Longitude       => Float64,
                :Elevation       => Float64,
                Symbol("d18O.Piso.AI")     => Float64,
                Symbol("d2H.Piso.AI")      => Float64)

    input_meteoiso = @linq DataFrame(File(path_meteoiso; header = 4,
                skipto=5, delim=',', ignorerepeated=true, types=parsing_types))
    select!(input_meteoiso, Not([:Column1]))

    expected_names = [String(k) for k in keys(parsing_types)]
    assert_colnames_as_expected(input_meteoiso, path_meteoiso, expected_names)

    #####
    # Special treatment of Piso.AI data: contains data on 15th of each month
    # -> a) only keep needed columns and rename them
    select!(input_meteoiso, Not([:Site, :Latitude, :Longitude, :Elevation]))
    rename!(input_meteoiso,
                Symbol("d18O.Piso.AI") => :d18O,
                Symbol("d2H.Piso.AI") => :d2H)
    # -> b) Shift all measured dates to end of month to be able to constnat interpolate backwards.
    #       Further repeat the very line concerning the first month in the timeseries to be
    #       able to interpolate between the first and last day of the first month.
    input_meteoiso_first = deepcopy(input_meteoiso[1,:])
    input_meteoiso_first.Date = floor(input_meteoiso_first.Date, Month)
    input_meteoiso = @linq input_meteoiso |>
            transform(Date = ceil.(:Date, Month))
    push!(input_meteoiso, input_meteoiso_first) # repeat the first
    sort!(input_meteoiso, [:Date])
    input_meteoiso[end,:].Date = input_meteoiso[end,:].Date - Day(1) # Modify last to remain within year (31.12. instead of 01.01)
    #####

    # Transform times from DateTimes to simulation time (Float of Days)
    input_meteoiso = @linq input_meteoiso |>
        transform(Date = DateTime2RelativeDaysFloat.(:Date, input_meteoveg_reference_date)) |>
        rename(Dict(:Date => :days))

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
        input_meteoiso = @linq input_meteoiso |>
            where(:days .>= startday, :days .<= endday)
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
function read_path_initial_conditions(path_initial_conditions)
    parsing_types =
        Dict(# Initial conditions (except for depth-dependent u_aux_PSIM) -------
            "u_GWAT_init_mm" => Float64,       "u_INTS_init_mm" => Float64,
            "u_INTR_init_mm" => Float64,       "u_SNOW_init_mm" => Float64,
            "u_CC_init_MJ_per_m2"   => Float64,       "u_SNOWLQ_init_mm" => Float64)
    input_initial_conditions = DataFrame(File(path_initial_conditions;
        transpose=true, drop=[1], comment = "###",
        types=parsing_types))

    expected_names = [String(k) for k in keys(parsing_types)]
    assert_colnames_as_expected(input_initial_conditions, path_initial_conditions, expected_names)

    return input_initial_conditions
end

# path_param = "examples/BEA2016-reset-FALSE-input/BEA2016-reset-FALSE_param.csv"
function read_path_param(path_param; simulate_isotopes::Bool = false)
    parsing_types =
        Dict(### Isotope tranpsport parameters  -------,NA
            "TODO" => Float64, "TODO2" => Float64,
            # Meteorologic site parameters -------
            "LAT_DEG" => Float64,
            "ESLOPE_DEG" => Float64,       "ASPECT_DEG" => Float64,
            "ALB" => Float64,              "ALBSN" => Float64,
            "C1" => Float64,               "C2" => Float64,               "C3" => Float64,
            "WNDRAT" => Float64,           "FETCH" => Float64,            "Z0W" => Float64,              "ZW" => Float64,
            # Canopy parameters -------
            "LWIDTH" => Float64,           "Z0G" => Float64,              "Z0S" => Float64,
            "LPC" => Float64,              "CS" => Float64,               "CZS" => Float64,
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

    input_param = DataFrame(File(path_param;
        transpose=true, drop=[1], comment = "###",
        types = parsing_types))

    expected_names = [String(k) for k in keys(parsing_types)]
    assert_colnames_as_expected(input_param, path_param, expected_names)

    # Set minimum/maximum values
    # from LWFBrook90R:PFILE.h
    input_param[:,:FXYLEM]  = min.(input_param[:,:FXYLEM], 0.990)
    input_param[:,:INITRLEN] = max.(input_param[:,:INITRLEN], 0.010)
    input_param[:,:INITRDEP] = max.(input_param[:,:INITRDEP], 0.010)

    return input_param
end

# path_storm_durations = "examples/BEA2016-reset-FALSE-input/BEA2016-reset-FALSE_meteo_storm_durations.csv"
function read_path_storm_durations(path_storm_durations)
    parsing_types =
        Dict("month" => String, "average_storm_duration_h" => Float64)
    input_storm_durations = DataFrame(File(path_storm_durations;
        strict=true, comment = "###",
        types=parsing_types))

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
        @error """
        Could not derive which hydraulic model parametrization (Mualem-Van-Genuchten or
        Clapp-Hornberger) to use, based the column names of the input file '$path_soil_horizons'.
        Please check and correct the input file!

        Expected column names are either:
        Mualem-Van-Genuchten: $MualVanGen_expected_column_names
        or for Clapp-Hornberger: $ClappHornberger_expected_column_names
        """
    end

    if FLAG_MualVanGen == 1
        # Assert units:
        assert_unitsHeader_as_expected(path_soil_horizons,
            DataFrame(HorizonNr = "-", Upper_m = "m", Lower_m = "m", ths_volFrac = "volume fraction (-)",
                thr_volFrac = "volume fraction (-)", alpha_perMeter = "perMeter", npar_ = "-",
                ksat_mmDay = "mm per day", tort_ = "-", gravel_volFrac = "volume fraction (-)"))

        # Load file
        input_soil_horizons = DataFrame(File(path_soil_horizons;
            datarow=3, header=MualVanGen_expected_column_names,
            types=[Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64],
            delim=','))
    else # FLAG_MualVanGen = 0 (Clapp-Hornberger)
        # Assert units:
        assert_unitsHeader_as_expected(path_soil_horizons,
            DataFram(HorizonNr = "-", Upper_m = "m", Lower_m = "m", thsat_volFrac = "volume fraction (-)",
                thetaf_volFrac = "volume fraction (-)", psif_kPa = "kPa", bexp_ = "-",
                kf_mmDay = "mm per day", wtinf_ = "-", gravel_volFrac = "volume fraction "))

        # Load file
        input_soil_horizons = DataFrame(File(path_soil_horizons;
            datarow=3, header=ClappHornberger_expected_column_names,
            types=[Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64],
            delim=','))# ignorerepeated=true
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

    return input_soil_horizons, FLAG_MualVanGen
end

# path_soil_discretization = "examples/BEA2016-reset-FALSE-input/BEA2016-reset-FALSE_soil_discretization.csv"
function read_path_soil_discretization(path_soil_discretization)
    parsing_types =
        Dict("Upper_m"      => Float64,
             "Lower_m"      => Float64,
             "Rootden_"      => Float64,
             "uAux_PSIM_init_kPa"   => Float64,
             "u_delta18O_init_mUr" => Float64,
             "u_delta2H_init_mUr"  => Float64)

    input_soil_discretization = DataFrame(File(path_soil_discretization;
        datarow=3, missingstring = "NA", types=parsing_types))

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
    # Assert units:
    assert_unitsHeader_as_expected(path_soil_discretization,
        DataFrame(Upper_m = "m", Lower_m = "m", Rootden_ = "-",
            uAux_PSIM_init_kPa = "kPa",
            u_delta18O_init_mUr = "mUr", u_delta2H_init_mUr = "mUr"))

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
    @assert DataFrame(File(path;datarow=2,limit=1)) == expected_units """
    Unexpected units in input file $path. Expected units:
    $expected_units"""
end

function assert_unitsColumn_as_expected(path, expected_units)
    # TODO(bernhard): define units in and implement this check.
    # @assert DataFrame(File(path;datarow=2,limit=1)) == expected_units """
    # Unexpected units in input file $path. Expected units:
    # $expected_units"""
end
