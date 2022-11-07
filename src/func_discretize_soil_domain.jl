# TODO: maybe move this into a module? on soil discretization?

"""
    discretize_soil(folder::String, prefix::String; suffix::String = "")

Load input file `soil_discretization.csv` for LWFBrook90:

The file `soil_discretization.csv` was created with an R script
`generate_LWFBrook90jl_Input.R` that takes the same arguements as the R function
`LWFBrook90R::run_LWFB90()` and generates the corresponding input files for
    LWFBrook90.jl.
"""
function discretize_soil(path_soil_discretization::String)
    # Load pre-defined discretization
    input_soil_discretization = LWFBrook90.read_path_soil_discretization(path_soil_discretization)

    return input_soil_discretization
end
"""
discretize_soil(;
    Δz_m::Vector{T},
    Rootden_::Function,
    uAux_PSIM_init_kPa::Function,
    u_delta18O_init_permil::Function = ((Δz_m) -> missing),
    u_delta2H_init_permil::Function = ((Δz_m) -> missing))

Manually generate a soil- and root-discretization for LWFBrook90.jl. This function can be
used as alternative to loading an input file `soil_discretization.csv`.

Reqired arguments are a vector with thickness of the discretization cells/layers in meter
`Δ_m`, and functions `Rootden_(Δ_m)` that generates the root density, a function
`uAux_PSIM_init_kPa(Δ_m)` that generates the initial values of ψ (kPa)
for all cells as a function of `Δ_m`.

Optinal arguments are functions initial conditions of the two isotopic
signatures based on `Δ_m`.
"""
function discretize_soil(;
    Δz_m::Vector{T},
    Rootden_::Function,
    uAux_PSIM_init_kPa::Function,
    u_delta18O_init_permil::Function = ((Δz_m) -> missing),
    u_delta2H_init_permil::Function = ((Δz_m) -> missing)) where {T<:Number}

    # Derive cell interfaces based on the cell spacing Δz_m
    z_cell_interfaces = -[0; cumsum(Δz_m)]
    z_cell_interfaces = round.(z_cell_interfaces; digits=6)

    # Output required DataFrame
    DataFrame(
        Upper_m = z_cell_interfaces[Not(end)],
        Lower_m = z_cell_interfaces[Not(1)],
        Rootden_ = Rootden_(Δz_m),
        uAux_PSIM_init_kPa = uAux_PSIM_init_kPa(Δz_m),
        u_delta18O_init_permil = u_delta18O_init_permil(Δz_m),
        u_delta2H_init_permil = u_delta2H_init_permil(Δz_m))
end

# function discretize(soil_horizons::DataFrame, )
#     # checks on structure
#     @assert names(soil_horizons) == ["HorizonNr", "Upper_m", "Lower_m", "shp"]
#     coltypes = [eltype(col) for col in eachcol(input_soil_horizons[1])]
#     @assert coltypes[1] <: Integer
#     @assert coltypes[2] <: Real
#     @assert coltypes[3] <: Real
#     @assert coltypes[4] <: LWFBrook90.AbstractSoilHydraulicParams

#     # return DataFrame(Upper_m, Lower_m, Rootden_, ICs, shp)
#     # where ICs is an Array{State}, and shp an Array{AbstractSoilHydraulicParams}
#     # TODO(bernhard): define
#     return DiscretizedSoilDomain(Upper_m, Lower_m, Rootden_, ICs, shp)
# end



# TODO(bernhard): unite this function with discretize_soil
"""
    refine_soil_discretization(
        input_soil_horizons,
        input_soil_discretization,
        soil_output_depths,
        IDEPTH_m,
        QDEPTH_m)

Discretize soil domain into computational layers and attribute soil parameters based on the defined horizons.
Densify discretization whenever an interface or an additional layer is needed. This is the case for interfaces
when a new soil horizons begins or an additional computational layer (consisting of upper and lower interfaces)
when a state variable needs to be extracted at a specified output depth.
"""
function refine_soil_discretization(
    input_soil_horizons,
    input_soil_discretization,
    soil_output_depths,
    IDEPTH_m,
    QDEPTH_m)

    # input_soil_horizons =  LWFBrook90.read_path_soil_horizons(
    #     "test/test-assets/DAV-2020/input-files/DAV_LW1_def_soil_horizons.csv");
    # input_soil_discretization = LWFBrook90.read_path_soil_discretization
    # This function maps the input parameters: (input_soil_horizons, ILAYER, QLAYER, INITRDEP, RGRORATE)
    # onto the soil discretization (input_soil_discretization)

    # input_soil_horizons: A matrix of the 8 soil materials parameters.
    # input_soil_discretization:
    # IDEPTH_m depth over which infiltration is distributed, [m] "IDEPTH_m determines the
    #        number of soil layers over which infiltration is distributed when INFEXP is
    #        greater than 0."
    # QDEPTH_m soil depth for SRFL calculation, 0 to prevent SRFL, [m] "QDEPTH_m determines the
    #        layers over which wetness is calculated to determine source area (SRFL)
    #        parameters."
    # INITRDEP
    # RGRORATE

    ############
    # 1) Check if discretization needs to be refined
    # 1a) Find out which layers need to be added
    layers_to_insert = soil_output_depths

    # ε = 0.001 # thickness of layer to be inserted, [m]
    # ε = 0.025 # thickness of layer to be inserted, [m]
    ε = 0.050 # thickness of layer to be inserted, [m]
    @assert all(abs.(diff(soil_output_depths)) .>= ε) """
        Requested soil_output_depths (additional layers) must be further away than $ε m.∇
        Requested were: $(soil_output_depths)
    """


    needed_interfaces_for_additional_layers = zeros(Float64, 0)
    for layer in layers_to_insert
        println(layer)
        if any(layer .∈ input_soil_discretization.Lower_m)
            # specific depth is already an interface in the discretization
            # -> only add a lower interface ε below the upper
            append!(needed_interfaces_for_additional_layers, round(layer - ε; digits=3))
        else
            # specific depth is not yet an interface
            # -> add upper and lower interface
            append!(needed_interfaces_for_additional_layers, round(layer;     digits=3))
            append!(needed_interfaces_for_additional_layers, round(layer - ε; digits=3))
            # TODO: would technically need to check that no interface between layer and layer + ε
        end
    end

    # 1b) Find out which interfaces need to be added
    all_needed_interfaces = unique(sort(
        [
         # interfaces for change in horizons or inflow depths or QDEPTH_m
         input_soil_horizons[!,"Upper_m"];
         input_soil_horizons[!,"Lower_m"];
         -IDEPTH_m; # m (positive) to m (negative)
         -QDEPTH_m;  # m (positive) to m (negative)

         # interfaces for layers determined in 1a)
         needed_interfaces_for_additional_layers
         ];
        rev = true))

    existing_interfaces = [0; input_soil_discretization[!,"Lower_m"]]
    to_add = all_needed_interfaces[(!).(all_needed_interfaces .∈ (existing_interfaces,))]
    if length(to_add) != 0
        @warn "Adding soil layers at depths $to_add, to allow for either requested IDEPTH/QDEPTH:  $IDEPTH_m/$QDEPTH_m, "*
              "or to accomodate defined soil horizons: $(unique(vcat(input_soil_horizons[!,"Upper_m"], input_soil_horizons[!,"Lower_m"])))"*
              ifelse(length(soil_output_depths)==0,".",", or `soil_output_depths`=$soil_output_depths.")
    end

    # Add them to the DataFrame
    soil_discretization = copy(input_soil_discretization) # otherwise input argument is modified in-place
    # soil_discretization = input_soil_discretization
    for interface_to_add in to_add
        for i in 1:nrow(soil_discretization)
            condition = interface_to_add > soil_discretization[i, "Lower_m"]
            #println("$i, interface_to_add:$interface_to_add, $(soil_discretization[i, "Lower_m"]) $condition")
            if (condition)
                interface_to_add
                pos = i+1
                insert!.(eachcol(soil_discretization),
                        pos,
                        Matrix(soil_discretization)[pos-1,:])
                soil_discretization[pos-1, "Lower_m"] = interface_to_add[1]
                soil_discretization[pos,   "Upper_m"] = interface_to_add[1]
                break
            elseif (i == nrow(soil_discretization))
                # If interface to add is below the lowest layer: add it only in case it is about ε lower
                if (interface_to_add >= soil_discretization[i, "Lower_m"] - ε)
                    pos = i+1
                    insert!.(eachcol(soil_discretization),
                        pos,
                        Matrix(soil_discretization)[pos-1,:])
                    soil_discretization[pos, "Upper_m"] = soil_discretization[pos-1, "Lower_m"]
                    soil_discretization[pos, "Lower_m"] = interface_to_add[1]
                    break
                end
            end
        end
    end

    # Define ILAYER and QLAYER (to be used internally instead of IDEPTH_m and QDEPTH_m)
    is_infiltration_layer_BOOLEAN = -IDEPTH_m .<= soil_discretization[!,"Lower_m"]
    is_SRFL_layer_BOOLEAN         = -QDEPTH_m .<= soil_discretization[!,"Lower_m"]
    ILAYER = sum(is_infiltration_layer_BOOLEAN) # lowest node where Bottom_m is below IDEPTH_m
    QLAYER = sum(is_SRFL_layer_BOOLEAN)         # lowest node where Bottom_m is below QDEPTH_m

    if (-IDEPTH_m < soil_discretization[end,"Lower_m"]) ||
        (-QDEPTH_m < soil_discretization[end,"Lower_m"])
        @error "QDEPTH_m or IDEPTH_m were defined deeper than the lowest simulation element."
    end
    ############

    ############
    # 2) Append soil horizon data to discretized soil domain
    # Assert expectations
    @assert input_soil_horizons[1,"Upper_m"] ≈ 0
    @assert soil_discretization[1,"Upper_m"] ≈ 0
    @assert input_soil_horizons[1,"Lower_m"] < 0
    @assert soil_discretization[1,"Lower_m"] < 0

    # Find for each soil_discretization node in which horizon it lies
   which_horizon = fill(0, nrow(soil_discretization))
    for i = 1:nrow(soil_discretization)
        # For each soil discretizations (i) check to which which horizon it belongs
        idx = findfirst(soil_discretization[i,"Lower_m"] .>= input_soil_horizons[:,"Lower_m"])

        if (isnothing(idx) && i == nrow(soil_discretization))
            idx = nrow(input_soil_horizons)
        end

        which_horizon[i] = input_soil_horizons[idx,"HorizonNr"]
    end
    soil_discretization[:,"HorizonNr"] = which_horizon

    # If lowest discretized layer is only ε lower than defined soil horizons extrapolate
    # the properties of the lowest soil horizon to that infinitesimal thin layer
    if ( (soil_discretization[end,"HorizonNr"] == 0) &
        (abs(soil_discretization[end,"Lower_m"] - input_soil_horizons[end,"Lower_m"]) < ε) )
        soil_discretization[end,"HorizonNr"] = nrow(input_soil_horizons)
    end

    nlayer_before_join = nrow(soil_discretization)
    soil_discretization = innerjoin(soil_discretization,
                                    input_soil_horizons, makeunique=true,
                                    #select(input_soil_horizons, Not([:Upper_m,:Lower_m])),
                                    on = :HorizonNr)
                                    # NOTE: not using leftjoin becaus it transforms types to Union{Missing, Float64} instead of only Float64
                                    # Therefore, use innerjoin and check manually that no rows in soil_discretization are lost
    @assert nrow(soil_discretization) == nlayer_before_join """
        When merging soil properties defined in input CSV-file containing soil horizons, some layers from
        the soil discretization were lost. This is most likely due to the soil horizons not covering the entire discretization domain.
        (Note that the discretized domain is defined in a separate input CSV-file.)

        Please check and correct the input files.
        """
        # (Note that the discretized domain is defined in a separate input CSV-file and is further extended to include internal variables IDEPTH_m and QDEPTH_m.)
    ############

    ############
    HEAT = 0 # flag for heat balance; not implemented; continuous_SPAC.params[:HEAT], @ hardcoded

    THICK_m        = soil_discretization[!,"Upper_m"] - soil_discretization[!,"Lower_m"] # thickness of soil layer [m]
    THICK          = 1000*(THICK_m)                                    # thickness of soil layer [mm]
    PSIM_init      = soil_discretization[!,"uAux_PSIM_init_kPa"]       # initial condition PSIM [kPa]
    d18O_soil_init = soil_discretization[!,"u_delta18O_init_permil"] # initial condition soil water δ18O [‰]
    d2H_soil_init  = soil_discretization[!,"u_delta2H_init_permil"]  # initial condition soil water δ2H [‰]

    @assert all(PSIM_init .<= 0) "Initial matrix psi must be negative or zero"

    NLAYER = nrow(soil_discretization)

    return Dict([
                ("NLAYER",NLAYER),
                ("ILAYER",ILAYER),
                ("QLAYER",QLAYER),
                ("THICK",THICK),
                ("PSIM_init",PSIM_init),
                ("d18O_init",d18O_soil_init),
                ("d2H_init",d2H_soil_init),
                ("SHP", soil_discretization.shp),
                ("final_Rootden_", soil_discretization[!,"Rootden_"])])
                # ("PAR",PAR),
                # ("STONEF",STONEF),
                #("HeatCapOld",HeatCapOld),
                #("TopInfT", TopInfT)])

    # BELOW IS ONLY UNUSED CODE:...
    for i = 1:NLAYER
        if (HEAT != 0)
            @error "HEAT must be set to zero, as heat balance is not implemented."
        end
        # if (HEAT == 1)
        #     # TemperatureNew(i)       = soil_discretization[i,7] we don't have it in the input file!!!
        #     if i > NLAYER
        #         MUE[i] = THICK[i] / ( THICK[i] + THICK(I+1) )
        #         ZL[i]  = 0.5 * ( THICK[i] + THICK(I+1) )
        #     else
        #         MUE[i] = 0.5
        #         ZL[i]  = THICK[i]
        #     end
        #     TMean[i]   = 0.
        # end
    end
    ############

    # ############
    # # 3) Bring to simple vector and matrix form

    # heat flow -------
    # nmat   = nrow(input_soil_horizons)
    if (HEAT != 0) @error "HEAT must be set to zero, as heat balance is not implemented." end
    #       if (HEAT == 1)
    #        READ (12,*) Comment
    #        READ (12,*) tTop, Comment
    #        tTop=tTop
    #        READ (12,*) tBot, Comment
    #        tBot=tBot
    #        READ (12,*) TopInfT, Comment
    #        READ (12,*) BotInfT, Comment
    #        READ (12,*) kTopT, Comment
    #        READ (12,*) kBotT, Comment
    #        DO 207 I = 1, 7
    #         READ (12,*) Comment
    # 207    CONTINUE
    #        DO 208 I = 1, nmat
    #         READ (12,*) ilay, SV[i], OV[i], HB1[i], HB2[i], HB3[i]
    #         TPar[1,I) = SV[i]
    #         TPar[2,I) = OV[i]
    #         TPar[3,I) = THDis
    # C        thermal conductivities -- transfer from [J m-1 s-1 K-1] to  [J mm-1 d-1 K-1]
    #         TPar[4,I) = HB1[i] * 86.400
    #         TPar[5,I) = HB2[i] * 86.400
    #         TPar[6,I) = HB3[i] * 86.400
    # C         volumetric heat capacities for solid, water and organic -- transfer from [MJ m-2 mm-1 K-1] to [J mm-3 K-1]
    #         TPar[7,I) = p_CVSOL   # To define in module_CONSTANTS.jl: p_CVSOL = ?  # CVSOL  - volumetric heat capacity of solid (MJ m-2 mm-1 K-1)
    #         TPar[8,I) = p_CVORG   # To define in module_CONSTANTS.jl: p_CVORG = ?  # volumetric heat capacity of organic material (MJ m-2 mm-1 K-1) (hillel98)
    #         TPar[9,I) = LWFBrook90.CONSTANTS.p_CVLQ
    # 208    CONTINUE
    #        READ (12,*) C
    #       end

    # mat       = input_soil_horizons[!,"mat"]
    # ### # from LWFBrook90R:md_brook90.f95 (lines 105)
    # TPar = fill(NaN, (nmat, 10))
    # TPar .= -1
    # HeatCapOld = fill(NaN, NLAYER)
    # for i = 1:NLAYER
    #     HeatCapOld[i] = TPar[mat[i],7] * TPar[mat[i],1] + TPar[mat[i],8]
    #                 # TPar[2,mat[i])+TPar[9,mat[i])*SWATI[i]/THICK[i]
    # end
    ###
end





"""
    Rootden_beta_(
    β;
    Δz_m,
    z_rootMax_m = maximum(cumsum(Δz_m)),
    z_Upper_m = 0)

Define the relative root density in each discretized soil layer based on the beta model from
(Gale and Grigal, 1987).

This function returns the instantaneous root fraction (dY/dd) (units of -/cm). Unless cropped
    by the total effective rooting depth the area under the curve sums up to 1.0 (when plotted vs cm).
    The method reduces the amount of roots in a discretization layer in case the effective maximal rooting
    depth comes to lie within that layer - it does not need to modify the discretization in that case.
"""
function Rootden_beta_(
    β;
    Δz_m,
    z_rootMax_m = -maximum(cumsum(Δz_m)), # total_effective_rooting_depth_m
    z_Upper_m = 0) # # Upper interface of topmost discretization cell
    # e.g. Δz_m = fill(0.10, 10)
    #      z_rootMax_m = -0.22
    #      β = 0.97


    # Notation:
    # z: absolute vertical position (positive upward, i.e. soil values are negative)
    # d: depth below beginning of mineral soil
    @assert z_Upper_m <= 0 "Discretization of humus layer (>0m) is untested."
    @assert minimum(Δz_m) > 0 "Discretization cell height Δz_m must be positive"
    z_interfaces_m = z_Upper_m .+ [0; cumsum(-Δz_m)]

    # (Gale and Grigal, 1987)
    Y_fct(β, d_cm) = 1 .- β .^ d_cm
    dYdd_fct(β, d_cm) = -β .^ d_cm .* log(β) # unused

    # # (Gale and Grigal, 1987) figures 1 and 2
    # d_midpoints_m = -z_Upper_m + Δz_m[1]/2 .+
    #     cumsum(Δz_m) .- Δz_m[1] # d is soil depth d in centimeters
    # # TODO: by using d_midpoint of we make a small error.
    # #       Better would be to use integral between cell boundaries.
    # d_cm = d_midpoints_m*100
    # pl_1 = plot(
    #     Y_fct(0.92, d_cm), d_cm,
    #     xflip = true,
    #     yflip = true, ylabel = "Depth_cm", xlabel = "Cumulative Root Fraction (Y)");
    # plot!(Y_fct(0.95, d_cm), d_cm);
    # plot!(Y_fct(0.97, d_cm), d_cm);
    # pl_2 = plot(
    #     dYdd_fct(0.92, d_cm), d_cm,
    #     yflip = true, ylabel = "Depth_cm", xlabel = "Instantaneous Root Fraction (dY/dd)");
    # plot!(dYdd_fct(0.95, d_cm), d_cm);
    # plot!(dYdd_fct(0.97, d_cm), d_cm);
    # plot(pl_1, pl_2; layout=(1,2))

    # 0) Preparatory steps
    # Depth of cell interfaces at which we evaluate the cumulative root fraction
    d_interfaces_cm = -z_interfaces_m * 100
    # # Crop cumulative density at total effective rooting depth
    d_root_depth_cm = -z_rootMax_m * 100

    # 1) Compute cumulative density inside of a cell (i.e. increase in cumulative denisity from upper to lower interface)
    rootden_increase_per_cell = diff(Y_fct(β, d_interfaces_cm))

    # 2) Correct for total effective rooting depth
    #   a) total rooting depth is below lowest discretization cell                    (nothing needs changing)
    #   b) total rooting depth coincides with lower boundary of a discretization cell (remove roots further below)
    #   c) total rooting depth is inside of a discretization cell                     (remove roots below and correct cell)
    d_upper_interfaces_cm = d_interfaces_cm[Not(end),]
    d_lower_interfaces_cm = d_interfaces_cm[Not(1),]

    # remove roots from cells that are completely below rooting depth
    rootden_increase_per_cell[d_upper_interfaces_cm.>=d_root_depth_cm] .= 0
    rootden_increase_per_cell

    # correct root density in cell where total rooting depth lies
    idx = findfirst(d_root_depth_cm .<= d_lower_interfaces_cm)
    if (!isnothing(idx))
        # First computation above was:
        # Y_fct(β, d_lower_interfaces_cm[idx]) - Y_fct(β, d_upper_interfaces_cm[idx])
        # Corrected computation is:
        rootden_increase_per_cell[idx] =
            Y_fct(β, d_root_depth_cm) - Y_fct(β, d_upper_interfaces_cm[idx])
        # (note if d_lower_interfaces_cm[idx] == d_root_depth_cm, this keeps the value before)
    end

    # 3) Normalize the increase in root density over the cell by the cell height
    #    note that this uses Δz_m also for the one cell which contains the total rooting depth
    rootden_increase_per_cmDepth = rootden_increase_per_cell ./ (Δz_m .* 100)
end

# Rootden_beta_(0.97, Δz_m = fill(0.10, 10))
# Rootden_beta_(0.97, Δz_m = fill(0.10, 10), z_rootMax_m = 0.2)

# Plot Rootden_beta_
# Δz_m = fill(0.10, 10) # (m)
# Δz_m = diff(collect(0:0.01:1)) # (m)
# z_rootMax_m = -0.225
# plot(Rootden_beta_(0.97, Δz_m = Δz_m),
#     cumsum(Δz_m * 100),
#     ylabel = "Depth (cm)", xlabel = "Instantaneous root fraction (dY/dd) (Area of 1.0, unless cropped)",
#     yflip = true, seriestype = :path, legend = :bottomright,
#     label = "No max Root")
# # plot!(Rootden_beta_(0.97, Δz_m = Δz_m, z_rootMax_m = z_rootMax_m),
# #     cumsum(Δz_m * 100), seriestype = [:scatter], label = "maxRoot: $(z_rootMax_m)m")
# # plot!(Rootden_beta_(0.97, Δz_m = Δz_m, z_rootMax_m = -1.2),
# #     cumsum(Δz_m * 100), seriestype = [:scatter], label = "maxRoot: $(-1.2)m")
# # plot!(Rootden_beta_(0.97, Δz_m = Δz_m, z_rootMax_m = -0.5),
# #     cumsum(Δz_m * 100), seriestype = [:scatter], label = "maxRoot: $(-0.5)m",)
# plot!(Rootden_beta_(0.97, Δz_m = Δz_m, z_rootMax_m = z_rootMax_m),
#     cumsum(Δz_m * 100), seriestype = [:path], label = "maxRoot: $(z_rootMax_m)m")
# plot!(Rootden_beta_(0.97, Δz_m = Δz_m, z_rootMax_m = -0.9),
#     cumsum(Δz_m * 100), seriestype = [:path], label = "maxRoot: $(-0.9)m")
# plot!(Rootden_beta_(0.97, Δz_m = Δz_m, z_rootMax_m = -0.5),
#     cumsum(Δz_m * 100), seriestype = [:path], label = "maxRoot: $(-0.5)m",)


"""
    HammelKennel_transient_root_density()

Take root density distribution parameters and generates root density in time (t) and space (z, negative downward).
"""
function HammelKennel_transient_root_density(;
        timepoints,   p_AGE,
        p_INITRDEP,   p_INITRLEN,
        p_RGROPER_y,  p_RGRORATE_m_per_y, # vertical root growth rate in m per y (accessing new soil layers)
        p_THICK,      final_Rootden_profile)

    NLAYER = length(p_THICK)
    # Hammel and Kennel 2001:

    # 0) preprocess final equilibrium state: (final_Rootden_profile)
    # final_Rootden_profile: relative values of final root density per unit volume
    # make root density continuous between first and last layer
    i1 = findfirst(final_Rootden_profile .> 1.e-6) # first layer where final_Rootden_profile[i] is >=1.e-6
    i2 = findlast(final_Rootden_profile .> 1.e-6)  # last layer (in 1:NLAYER) where final_Rootden_profile[i] is >=1.e-6
    if !isnothing(i1) && isnothing(i2)
        i2 = NLAYER
    end
    for i = i1:i2 # set final_Rootden_profile to minimum of 1.e-6 within the root zone (between i1 and i2)
        final_Rootden_profile[i] = max( final_Rootden_profile[i], 1.01e-6)
    end

    # 1) Define age at which root growth arrives in a specific depth and lateral root growth commences
    # z_center     = soil_discretization[!,"Upper_m"] - soil_discr["THICK"]/1000/2 # soil depth [m] (midpoint, i.e. average depth of layer)
    z_center_m     = -(cumsum(p_THICK) - p_THICK/2)/1000
    depfirst_m  = z_center_m[1] - p_THICK[1] / 1000.
    # Here z_center_m and depfirst_m are negative values in the soil. (z axis upward)

    # 1a) compute arrival times tini
    vᵣ = p_RGRORATE_m_per_y # growth rate (m/y)
    zᵣ₀ = p_INITRDEP # initial depth of root
    # zᵣ =  # maximal depth of root
    # root_arrival_time = (zᵣ - zᵣ₀)/vᵣ
    # tᵣ =  # time period of lateral growth (after first vertical arrival)
    # rₗ₀ = #
    # # t_ini(z) : arrival time of roots at depth z, then starts lateral growth

    t_ini = fill(NaN, NLAYER)
    for i = 1:NLAYER
        t_ini[i] = 1.e+20               # start out with "never"
        z = depfirst_m - z_center_m[i]
        roots_present = final_Rootden_profile[i] >= 1.e-6
        if roots_present
            if z <= zᵣ₀
                t_ini[i] = 0.
            elseif z > zᵣ₀
                if vᵣ > 0
                    t_ini[i] = (z - zᵣ₀)/vᵣ
                end
            end
        end
    end

    # 2) lateral growth per layer starts at t_ini and ends at t_ini + tᵣ
    # Time dependent root density parameters is a vector quantity that is dependent on time (i.e. 2D array):
    p_RELDEN_2Darray = fill(NaN, length(timepoints), NLAYER)
    for i in eachindex(timepoints)
        p_RELDEN_2Darray[i,:] =
            HammelKennel_lateral_rootgrowth(;
                final_Rootden_profile = final_Rootden_profile,
                INITRDEP_m        = p_INITRDEP,
                INITRLEN_m_per_m2 = p_INITRLEN,
                t_y               = p_AGE(timepoints[i]),
                tstart_y          = t_ini,
                RGROPER_yrs       = p_RGROPER_y)
    end
    p_fT_RELDEN =  extrapolate(interpolate((timepoints, 1:NLAYER), p_RELDEN_2Darray,
                                        (Gridded(Constant{Next}()), NoInterp()), # 1st dimension: ..., 2nd dimension NoInterp()
                                        ), Flat()) # extrapolate flat, alternative: Throw()
    return p_fT_RELDEN
end