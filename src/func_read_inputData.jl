using CSV: read, File
using DataFrames: DataFrame, rename# ,select
using DataFramesMeta#: @linq, transform, DataFramesMeta
using Dates: DateTime, Millisecond, Second, Day, month, value, dayofyear

"""
    read_LWFBrook90R_inputData(folder::String, prefix::String)

Load different input files for LWFBrook90:
- meteoveg.csv
- param.csv
- pdur.csv
- initial_conditions.csv
- soil_horizons.csv
- soil_discretization.csv

These files were created with an R script `generate_LWFBrook90Julia_Input.R` that
takes the same arguements as the R funciton `LWFBrook90R::run_LWFB90()` and generates
the corresponding Julia input functions.
"""
function read_LWFBrook90R_inputData(folder::String, prefix::String)
    # https://towardsdatascience.com/read-csv-to-data-frame-in-julia-programming-lang-77f3d0081c14
    ## A) Define paths of all input files
    path_meteoveg           = joinpath(folder, prefix*"_meteoveg.csv")
    path_param              = joinpath(folder, prefix*"_param.csv")
    # unused path_precdat=joinpath(folder, prefix*"_precdat.csv")
    path_pdur               = joinpath(folder, prefix*"_pdur.csv")
    path_initial_conditions = joinpath(folder, prefix*"_initial_conditions.csv")
    path_soil_horizons      = joinpath(folder, prefix*"_soil_horizons.csv")
    path_soil_discretization= joinpath(folder, prefix*"_soil_discretization.csv")

    ## B) Load input data (time- and/or space-varying parameters)
    input_meteoveg, input_meteoveg_reference_date = read_path_meteoveg(path_meteoveg)

    ## C) Load other parameters
    # Load model input parameters
    #' @param param A numeric vector of model input parameters. Order:
    input_param = read_path_param(path_param)

    # Load precipitation data
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

    #' @param pdur a [1,12]-matrix of precipitation durations (hours) for each month.
    input_pdur = read_path_pdur(path_pdur)

    # Load initial conditions of scalar state variables
    input_initial_conditions = read_path_initial_conditions(path_initial_conditions)

    # Load soil data

    #' @param soil_materials A matrix of the 8 soil materials parameters.
    input_soil_horizons = read_path_soil_horizons(
        path_soil_horizons,
        input_param[1,"FLAG_MualVanGen"])
    #' @param soil_nodes A matrix of the soil model layers with columns
    #'                   nl (layer number), layer midpoint (m), thickness (mm), mat, psiini (kPa), rootden (-).
    input_soil_discretization = read_path_soil_discretization(path_soil_discretization)

    return (input_meteoveg,
            input_meteoveg_reference_date,
            input_param,
            input_pdur,
            input_initial_conditions,
            input_soil_horizons,
            input_soil_discretization)
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

function read_path_meteoveg(path_meteoveg)
    # Load meteo
    input_meteoveg = @linq DataFrame(File(path_meteoveg;
        datarow=2, delim=',', ignorerepeated=true,
        header=["dates", #"YY","MM","DD",
                "GLOBRAD","TMAX","TMIN","VAPPRES","WIND",
                "PRECIN","DENSEF","HEIGHT","LAI","SAI","AGE"],
        types=Dict(:dates => DateTime, #YY,MM,DD,
                :GLOBRAD => Float64,:TMAX => Float64,:TMIN => Float64,:VAPPRES => Float64,:WIND => Float64,
                :PRECIN => Float64,:DENSEF => Float64,:HEIGHT => Float64,:LAI => Float64,:SAI => Float64,:AGE => Float64))) |>
        transform(dates = DateTime.(:dates))

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

function read_path_initial_conditions(path_initial_conditions)
    input_initial_conditions = DataFrame(File(path_initial_conditions;
        transpose=true, drop=[1], comment = "###",
        types = Dict(# Initial conditions (except for depth-dependent u_aux_PSIM) -------
                    "u_GWAT_init" => Float64,       "u_INTS_init" => Float64,
                    "u_INTR_init" => Float64,       "u_SNOW_init" => Float64,
                    "u_CC_init" => Float64,         "u_SNOWLQ_init" => Float64)))
    # Assert that inputs are correct:
    received_names = names(input_initial_conditions)
    expected_names = [
        "u_GWAT_init","u_INTS_init","u_INTR_init","u_SNOW_init","u_CC_init","u_SNOWLQ_init",
        ]

    @assert all(expected_names .∈ (received_names,)) """
    Input file '$path_initial_conditions' does not contain all the expected entries.
    \nExpected but missing:\n$(expected_names[(!).(expected_names .∈ (received_names,))])
    \n\nExpected:\n$expected_names\nReceived:\n$received_names\n
    """
    @assert all(received_names .∈ (expected_names,)) """
    Input file '$path_initial_conditions' contains other than the expected entries.
    \nReceived unexpected: $(received_names[(!).(received_names .∈ (expected_names,))])
    \nExpected:\n$expected_names\nReceived:\n$received_names\n
    """

    return input_initial_conditions
end

function read_path_param(path_param)
    input_param = DataFrame(File(path_param;
        transpose=true, drop=[1], comment = "###",
        types = Dict(# Meteorologic site parameters -------
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
                    "MAXRTLN" => Float64,          "INITRLEN" => Float64,         "INITRDEP" => Float64,
                    "RGRORATE" => Float64,         "RGROPER" => Float64,          "FXYLEM" => Float64,
                    "PSICR" => Float64,            "RTRAD" => Float64,            "NOOUTF" => Int64,    # TODO(bernhard): make boolean
                    # Soil parameters -------
                    "IDEPTH" => Int64,             "QDEPTH" => Int64,
                    "FLAG_MualVanGen" => Int64,    # TODO(bernhard): make boolean
                    "RSSA" => Float64,             "RSSB" => Float64,             "INFEXP" => Float64,
                    "BYPAR" => Int64  ,             # TODO(bernhard): make boolean
                    "QFPAR" => Float64,            "QFFC" => Float64,
                    "IMPERV" => Float64,           "DSLOPE" => Float64,           "LENGTH_SLOPE" => Float64,
                    "DRAIN" => Float64,            "GSC" => Float64,              "GSP" => Float64,
                    # Numerical solver parameters -------
                    "DTIMAX" => Float64,           "DSWMAX" => Float64,           "DPSIMAX" => Float64)))

    # Assert that inputs are correct:
    received_names = names(input_param)
    expected_names = [
        "LAT_DEG","ESLOPE_DEG","ASPECT_DEG",
        "ALB","ALBSN","C1","C2","C3","WNDRAT","FETCH","Z0W","ZW",
        "LWIDTH","Z0G","Z0S","LPC","CS","CZS","CZR","HS","HR","ZMINH","RHOTP","NN",
        "FRINTLAI","FSINTLAI","FRINTSAI","FSINTSAI","CINTRL",
        "CINTRS","CINTSL","CINTSS","RSTEMP","MELFAC","CCFAC","LAIMLT","SAIMLT","GRDMLT",
        "MAXLQF","KSNVP","SNODEN","GLMAX","CR","GLMIN","RM","R5","CVPD","TL","T1","T2",
        "TH","MXKPL","MXRTLN","INITRLEN","INITRDEP","RGRORATE","RGROPER","FXYLEM","PSICR",
        "RTRAD","NOOUTF","IDEPTH","QDEPTH","FLAG_MualVanGen","RSSA","RSSB","INFEXP","BYPAR",
         "QFPAR","QFFC","IMPERV","DSLOPE","LENGTH_SLOPE","DRAIN","GSC","GSP",
         "DTIMAX","DSWMAX","DPSIMAX"]

    @assert all(expected_names .∈ (received_names,)) """
    Input file '$path_param' does not contain all the expected entries.
    \nExpected but missing:\n$(expected_names[(!).(expected_names .∈ (received_names,))])
    \n\nExpected:\n$expected_names\nReceived:\n$received_names\n
    """
    @assert all(received_names .∈ (expected_names,)) """
    Input file '$path_param' contains other than the expected entries.
    \nReceived unexpected: $(received_names[(!).(received_names .∈ (expected_names,))])
    \nExpected:\n$expected_names\nReceived:\n$received_names\n
    """

    # from LWFBrook90R:PFILE.h
    input_param[:,:FXYLEM]  = min.(input_param[:,:FXYLEM], 0.990)
    input_param[:,:INITRLEN] = max.(input_param[:,:INITRLEN], 0.010)
    input_param[:,:INITRDEP] = max.(input_param[:,:INITRDEP], 0.010)

    return input_param
end

function read_path_pdur(path_pdur)
    input_pdur = DataFrame(File(path_pdur;
                           types=[String, Int64], strict=true, comment = "###",
                           datarow=2, header=["month", "pdur_h"], delim=',')) # ignorerepeated=true
    received_month_names = input_pdur[!,"month"]
    expected_month_names = ["January", "Februrary", "March", "April", "May", "June", "July",
                            "August", "September", "October", "November", "December"]
    @assert received_month_names == expected_month_names """
        Input file '$path_pdur' does not contain the months in the expected order.
        Please correct the input file.
        \n\nExpected:\n$expected_month_names\nReceived:\n$received_month_names\n
        """

    return input_pdur
end

function read_path_soil_horizons(path_soil_horizons, FLAG_MualVanGen)
    if FLAG_MualVanGen == 1
        input_soil_horizons = DataFrame(File(path_soil_horizons;
            # datarow=2, header=["HorizonNr","Upper_m","Lower_m","ths","thr","alpha","npar","ksat","tort","gravel"],
            datarow=2, header=["HorizonNr","Upper_m","Lower_m","ths_volFrac","thr_volFrac","alpha_perMeter","npar_","ksat_mmDay","tort_","gravel_volFrac"],
            types            =[Int64,      Float64,  Float64,  Float64,      Float64,      Float64,         Float64, Float64,    Float64,Float64],
            delim=','))# ignorerepeated=true
    else # FLAG_MualVanGen = 0 (Clapp-Hornberger)
        input_soil_horizons = DataFrame(File(path_soil_horizons;
            datarow=2, header=["HorizonNr","Upper_m","Lower_m","ths_volFrac","thr_volFrac","alpha_perMeter","npar_","ksat_mmDay","tort_","gravel_volFrac"],
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

    return input_soil_horizons
end

function read_path_soil_discretization(path_soil_discretization)
    input_soil_discretization = DataFrame(File(path_soil_discretization;
                                 datarow=2, header=["Upper_m","Lower_m","Rootden","Psiini_kPa","delta18O_mUr","delta2H_mUr"],
                                 types =           [Float64,  Float64,  Float64,  Float64,      Float64,      Float64],
                                 delim=','))# ignorerepeated=true

    # Check that defined layers do not overlap
    @assert input_soil_discretization[1:end-1,"Lower_m"] == input_soil_discretization[2:end,"Upper_m"] """
        Input file '$path_soil_discretization' contains overlapping layers.
        Please check and correct the input file.
        """

    @assert input_soil_discretization[1,"Upper_m"] ≈ 0 """
        Input file '$path_soil_discretization' does not start with 0.0 as Upper_m limit of first layer.
        Please check and correct the input file.
        """
    return input_soil_discretization
end