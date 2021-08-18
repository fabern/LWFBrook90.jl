using CSV: read, File
using DataFrames: DataFrame, rename# ,select
using DataFramesMeta#: @linq, transform, DataFramesMeta
using Dates: DateTime, Millisecond, Second, Day, month, value, dayofyear

"""
    read_LWFBrook90R_inputData(folder::String, prefix::String)

Load different input files for LWFBrook90:
- meteoveg
- param
- siteparam
- pdur
- soil_materials.csv
- soil_nodes.csv

These files were created with an R script `generate_LWFBrook90Julia_Input.R` that
takes the same arguements as the R funciton `LWFBrook90R::run_LWFB90()` and generates
the corresponding Julia input functions.
"""
function read_LWFBrook90R_inputData(folder::String, prefix::String)
    # https://towardsdatascience.com/read-csv-to-data-frame-in-julia-programming-lang-77f3d0081c14

    ## A) Define paths of all input files
    # TODO(bernhard): prepend path_
    path_meteoveg       = joinpath(folder, prefix*"_meteoveg.csv")
    path_param          = joinpath(folder, prefix*"_param.csv")
    path_siteparam      = joinpath(folder, prefix*"_siteparam.csv")
    # unused path_precdat=joinpath(folder, prefix*"_precdat.csv")
    path_pdur           = joinpath(folder, prefix*"_pdur.csv")
    path_soil_materials = joinpath(folder, prefix*"_soil_materials.csv")
    path_soil_nodes     = joinpath(folder, prefix*"_soil_nodes.csv")

    ## B) Load input data (time- and/or space-varying parameters)
    input_meteoveg, input_meteoveg_reference_date = read_path_meteoveg(path_meteoveg)

    ## C) Load other parameters
    # Load model input parameters
    #' @param param A numeric vector of model input parameters. Order:
    input_param = read_path_param(path_param)

    # Load site parameters
    #' @param siteparam A [1,6] matrix with site level information:
    #'                  start year, start doy, latitude, initial condition snow, initial condition groundwater, precipitation interval.
    input_siteparam = read_path_siteparam(path_siteparam)

    # Load precipitation data
    #' @param precdat A matrix of precipitation interval data with 6 columns:
    #'                   year, month, day, interval-number (1:precint), prec, mesflp.
    input_precdat = DataFrame(a = Nothing, b = Nothing)
    # TODO(bernhard): currently not implemented.
    #                 Only using PRECIN (resolution daily).
    #                 PRECDAT would allow to have smaller resolution (would require changes).
    # unused: input_precdat = read(path_precdat, DataFrame;
    # unused:           missingstring = "-999",
    # unused:           datarow=2, header=["year","month","day","interval_number","prec","mesflp"], delim=',',
    # unused:           ignorerepeated=true)

    #' @param pdur a [1,12]-matrix of precipitation durations (hours) for each month.
    input_pdur = read_path_pdur(path_pdur)

    # Load soil data

    #' @param soil_materials A matrix of the 8 soil materials parameters.
    #'                       When FLAG_MualVanGen = 1 (Mualem-van Genuchten), these refer to:
    #'                             mat, ths, thr, alpha (m-1), npar, ksat (mm d-1), tort (-), stonef (-).
    #'                       When FLAG_MualVanGen = 0 (Clapp-Hornberger):
    #'                             mat, thsat, thetaf, psif (kPa), bexp, kf (mm d-1), wetinf (-), stonef (-).
    input_soil_materials = DataFrame(File(path_soil_materials;
                                     types=[Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64],
                                     datarow=2, header=["mat","ths","thr","alpha","npar","ksat","tort","gravel"],
                                     delim=','))# ignorerepeated=true

    #' @param soil_nodes A matrix of the soil model layers with columns
    #'                   nl (layer number), layer midpoint (m), thickness (mm), mat, psiini (kPa), rootden (-).
    input_soil_nodes = DataFrame(File(path_soil_nodes;
                                 types=[Int64, Float64, Float64, Int64, Float64, Float64],
                                 datarow=2, header=["layer","midpoint","thick","mat","psiini","rootden"],
                                 delim=','))# ignorerepeated=true
    return (input_meteoveg,
            input_meteoveg_reference_date,
            input_param,
            input_siteparam,
            input_precdat,
            input_pdur,
            input_soil_materials,
            input_soil_nodes)
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
                "PRECIN","MESFL","DENSEF","HEIGHT","LAI","SAI","AGE"],
        types=Dict(:dates => DateTime, #YY,MM,DD,
                :GLOBRAD => Float64,:TMAX => Float64,:TMIN => Float64,:VAPPRES => Float64,:WIND => Float64,
                :PRECIN => Float64,:MESFL => Float64,:DENSEF => Float64,:HEIGHT => Float64,:LAI => Float64,:SAI => Float64,:AGE => Float64))) |>
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

function read_path_param(path_param)
    # input_param = DataFrame(File(path_param; types=[String, Float64], strict=true))
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
                    # Interception initial conditions -------
                    "u_INTR_init" => Float64,       "u_INTS_init" => Float64,
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
                    "PSICR" => Float64,            "RTRAD" => Float64,            "NOOUTF" => Int64,    # TODO(bernhard): make boolena
                    # Soil parameters -------
                    "ILAYER" => Int64  ,            # TODO(bernhard): switch to depth in mm
                    "QLAYER" => Int64  ,            # TODO(bernhard): switch to depth in mm
                    "FLAG_MualVanGen" => Int64  ,            # TODO(bernhard): make boolena
                    "RSSA" => Float64,             "RSSB" => Float64,             "INFEXP" => Float64,
                    "BYPAR" => Int64  ,             # TODO(bernhard): make boolena
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
        "u_INTR_init","u_INTS_init","FRINTLAI","FSINTLAI","FRINTSAI","FSINTSAI","CINTRL",
        "CINTRS","CINTSL","CINTSS","RSTEMP","MELFAC","CCFAC","LAIMLT","SAIMLT","GRDMLT",
        "MAXLQF","KSNVP","SNODEN","GLMAX","CR","GLMIN","RM","R5","CVPD","TL","T1","T2",
        "TH","MXKPL","MXRTLN","INITRLEN","INITRDEP","RGRORATE","RGROPER","FXYLEM","PSICR",
        "RTRAD","NOOUTF","ILAYER","QLAYER","FLAG_MualVanGen","RSSA","RSSB","INFEXP","BYPAR",
         "QFPAR","QFFC","IMPERV","DSLOPE","LENGTH_SLOPE","DRAIN","GSC","GSP",
         "DTIMAX","DSWMAX","DPSIMAX"]

    @assert all(expected_names .∈ (received_names,)) """
    Input file '$path_param' does not contain all the expected 80 entries.
    \nExpected but missing:\n$(expected_names[(!).(expected_names .∈ (received_names,))])
    \n\nExpected:\n$expected_names\nReceived:\n$received_names\n
    """
    @assert all(received_names .∈ (expected_names,)) """
    Input file '$path_param' contains other than the expected 80 entries.
    \nReceived unexpected: $(received_names[(!).(received_names .∈ (expected_names,))])
    \nExpected:\n$expected_names\nReceived:\n$received_names\n
    """

    # from LWFBrook90R:PFILE.h
    input_param[:,:FXYLEM]  = min.(input_param[:,:FXYLEM], 0.990)
    input_param[:,:INITRLEN] = max.(input_param[:,:INITRLEN], 0.010)
    input_param[:,:INITRDEP] = max.(input_param[:,:INITRDEP], 0.010)

    return input_param
end

function read_path_siteparam(path_siteparam)
    DataFrame(File(path_siteparam;
                            ##datarow=2, header=["start_year","start_doy","precip_interval_NPINT","LAT_DEG","u_SNOW_init","u_GWAT_init"],
                            #types=[Int64, Int64, Int64, Float64, Float64, Float64], strict=true,
                            datarow=2, header=["precip_interval_NPINT","u_SNOW_init","u_GWAT_init"],
                            types=[Int64, Float64, Float64], strict=true,
                            delim=',')) # ignorerepeated=true
end
function read_path_pdur(path_pdur)
    input_pdur = DataFrame(File(path_pdur;
                           types=[String, Int64], strict=true,
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