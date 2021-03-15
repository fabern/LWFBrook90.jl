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
    # Load meteo
    input_meteoveg = @linq DataFrame(File(path_meteoveg;
                            datarow=2, delim=',', ignorerepeated=true,
                            header=["dates", #"YY","MM","DD",
                                    "GLOBRAD","TMAX","TMIN","VAPPRES","WIND",
                                    "PRECIN","MESFL","DENSEF","HEIGHT","LAI","SAI","AGE"])) |>
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

    ## C) Load other parameters
    # Load model input parameters
    #' @param param A numeric vector of model input parameters. Order (derived from param_to_rlwfbrook90()):
    input_param = DataFrame(File(path_param; types=[String, Float64], strict=true))

    # Load site parameters
    #' @param siteparam A [1,6] matrix with site level information:
    #'                  start year, start doy, latitude, initial condition snow, initial condition groundwater, precipitation interval.
    input_siteparam = DataFrame(File(path_siteparam;
                                types=[Int64, Int64, Float64, Float64, Float64, Int64], strict=true,
                                datarow=2, header=["start_year","start_doy","lat","SNOW_init","GWAT_init","precip_interval"],
                                delim=',')) # ignorerepeated=true
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
    input_pdur = DataFrame(File(path_pdur;
                           types=[Int64], strict=true,
                           datarow=2, header=["pdur_h"], delim=',')) # ignorerepeated=true

    # Load soil data

    #' @param soil_materials A matrix of the 8 soil materials parameters.
    #'                       When imodel = 1 (Mualem-van Genuchten), these refer to:
    #'                             mat, ths, thr, alpha (m-1), npar, ksat (mm d-1), tort (-), stonef (-).
    #'                       When imodel = 2 (Clapp-Hornberger):
    #'                             mat, thsat, thetaf, psif (kPa), bexp, kf (mm d-1), wetinf (-), stonef (-).
    input_soil_materials = DataFrame(File(path_soil_materials;
                                     types=[Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64],
                                     datarow=2, header=["mat","ths","thr","alpha","npar","ksat","tort","gravel"],
                                     delim=','))# ignorerepeated=true

    #' @param soil_nodes A matrix of the soil model layers with columns
    #'                   nl (layer number), layer midpoint (m), thickness (mm), mat, psiini (kPa), rootden (-).
    input_soil_nodes = DataFrame(File(path_soil_nodes;
                                 types=[Int64,Float64, Float64, Int64, Float64, Float64],
                                 datarow=2, header=["layer","midpoint","thick","mat","psiini","rootden"],
                                 delim=','))# ignorerepeated=true

    return (input_meteoveg,
            input_param,
            input_siteparam,
            input_precdat,
            input_pdur,
            input_soil_materials,
            input_soil_nodes,
            reference_date)
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