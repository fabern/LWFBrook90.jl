############################################################################################
############################################################################################
############################################################################################
# Helpers for integration testing

function RMS_differences(sim, ref)
    # computes root mean squared differences

    if ("time" in names(sim) && "time" in names(ref))
        # case where: typeof(sim) <: DataFrame
        differences = Matrix(sim[:,Not(:time)]) .- Matrix(ref[:,Not(:time)])
    else # case where: typeof(sim) <: Matrix
        differences = sim .- ref
    end
    mean_squared = sum(differences.^2) / length(differences)

    return sqrt(mean_squared)
end

function prepare_θψδ_from_sim_and_reference(;
    path_jl_prefix, path_R_layeroutput, path_Hydrus)
    # path_jl_prefix      = "test-assets/Hammel-2001/input-files/Hammel_loam-NLayer-27-RESET=FALSE"
    # path_R_layeroutput  = "test-assets/Hammel-2001/output_LWFBrook90R/Hammel_loam-NLayer-27-RESET=TRUE_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv"
    # path_Hydrus         = "test-assets/Hammel-2001/output_Hydrus1D/Hammel_Test_Loam"

    # Run  simulation
    sim_sol, _, _ = run_simulation([dirname(path_jl_prefix), basename(path_jl_prefix)]) # "Hammel_loam-NLayer-27-RESET=TRUE"]

    # Load reference solution
    ## Load LWFbrook90R solution
    # referenceSolution_daily = DataFrame(File(path_R_dailyoutput))
    referenceSolution_layer = DataFrame(File(path_R_layeroutput))

    ## Load Hydrus solution
    HydrusSolution_denseTime = DataFrame(File(joinpath(path_Hydrus,"Obs_Node_processed.csv")))
    HydrusSolution_sparseTime = DataFrame(File(joinpath(path_Hydrus,"Nod_inf_processed.csv");
         missingstrings = ["","NA"]
                ))[:,[:timeStep, :depth_cm, :theta, :head_cm, :delta18O, :delta2H]]
    HydrusSolution_sparseTime = rename(subset(HydrusSolution_sparseTime,
        :depth_cm => ByRow(in([-10., -50., -100., -150., -190.]))),
        :timeStep => :time)
    # correct units of head:
    cm_to_kPa = 1/100*9.81*998/1000
    HydrusSolution_sparseTime.head_kPa = HydrusSolution_sparseTime.head_cm*cm_to_kPa

    # Compare simulation and reference
    ## Extract various variables at certain depths
    depth_to_read_out_mm = [100 500 1000 1500 1900]
    idx = find_indices(depth_to_read_out_mm, sim_sol)

    ### Sim
    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
            get_auxiliary_variables(sim_sol)
    # (u_δ18O_soil, u_δ2H_soil) =
    #         get_δ_soil(sim_sol)
    sim_θ = DataFrame(u_aux_θ[:,idx], :auto)
    sim_θ.time = sim_sol.t

    sim_ψ = DataFrame(u_aux_PSIM[:,idx], :auto)
    sim_ψ.time = sim_sol.t

    # sim_δ18O = DataFrame(u_δ18O_soil[:,idx], :auto)
    # sim_δ18O.time = sim_sol.t
    sim_δ18O = allowmissing(copy(sim_ψ))
    sim_δ18O[:,Not(:time)] .= missing

    # sim_δ2H = DataFrame(u_δ2H_soil[:,idx], :auto)
    # sim_δ2H.time = sim_sol.t
    sim_δ2H = allowmissing(copy(sim_ψ))
    sim_δ2H[:,Not(:time)] .= missing

    sim = (θ = sim_θ, ψ = sim_ψ, δ18O = sim_δ18O, δ2H = sim_δ2H)

    ### Ref
    ref_θ = unstack(
        referenceSolution_layer[:,[:yr, :mo, :da, :doy, :nl, :theta]],
        [:yr, :mo, :da, :doy], # ID variable
        :nl,   # variable that will be spread out into column names
        :theta,# value that will be filled into the wide table
        renamecols=x->Symbol(:nl_, x))[:, [:doy; Symbol.("nl_" .* string.(idx))]]
    rename!(ref_θ, :doy => :time)
    ref_ψ = unstack(
        referenceSolution_layer[:,[:yr, :mo, :da, :doy, :nl, :psimi]],
        [:yr, :mo, :da, :doy], # ID variable
        :nl,   # variable that will be spread out into column names
        :psimi,# value that will be filled into the wide table
        renamecols=x->Symbol(:nl_, x))[:, [:doy; Symbol.("nl_" .* string.(idx))]]
    rename!(ref_ψ, :doy => :time)

    ref_δ18O = allowmissing(copy(ref_ψ))
    ref_δ18O[:,Not(:time)] .= missing

    ref_δ2H = allowmissing(copy(ref_ψ))
    ref_δ2H[:,Not(:time)] .= missing

    ### Hydrus
    HydrusSolution_denseTime_θ_wide = unstack(
        HydrusSolution_denseTime[:, [:time, :depth, :theta]],
        [:time], # ID variable
        :depth,  # variable that will be spread out into column names
        :theta,  # value that will be filled into the wide table
        renamecols=x->Symbol(:depth_, x)
        )
    HydrusSolution_sparseTime_θ_wide = unstack(
        HydrusSolution_sparseTime[:, [:time, :depth_cm, :theta]],
        [:time],    # ID variable
        :depth_cm,  # variable that will be spread out into column names
        :theta,     # value that will be filled into the wide table
        renamecols=x->Symbol(:depth_cm_, x)
        )
    HydrusSolution_sparseTime_ψ_wide = unstack(
        HydrusSolution_sparseTime[:, [:time, :depth_cm, :head_kPa]],
        [:time],    # ID variable
        :depth_cm,  # variable that will be spread out into column names
        :head_kPa,   # value that will be filled into the wide table
        renamecols=x->Symbol(:depth_cm_, x)
        )
    HydrusSolution_sparseTime_δ2H_wide = unstack(
        HydrusSolution_sparseTime[:, [:time, :depth_cm, :delta2H]],
        [:time],  # ID variable
        :depth_cm,    # variable that will be spread out into column names
        :delta2H,     # value that will be filled into the wide table
        renamecols=x->Symbol(:depth_cm_, x)
        )
    HydrusSolution_sparseTime_δ18O_wide = unstack(
        HydrusSolution_sparseTime[:, [:time, :depth_cm, :delta18O]],
        [:time],  # ID variable
        :depth_cm,    # variable that will be spread out into column names
        :delta18O,     # value that will be filled into the wide table
        renamecols=x->Symbol(:depth_cm_, x)
        )

    # return θ from hydrus in daily resolution
    hyd_θ = HydrusSolution_sparseTime_θ_wide[HydrusSolution_sparseTime_θ_wide[:, :time] .% 1.0 .== 0, :] # Not(:time)]
    hyd_ψ = HydrusSolution_sparseTime_ψ_wide[HydrusSolution_sparseTime_ψ_wide[:, :time] .% 1.0 .== 0, :] # Not(:time)]
    hyd_δ18O = HydrusSolution_sparseTime_δ18O_wide[HydrusSolution_sparseTime_δ18O_wide[:, :time] .% 1.0 .== 0, :] # Not(:time)]
    hyd_δ2H = HydrusSolution_sparseTime_δ2H_wide[HydrusSolution_sparseTime_δ2H_wide[:, :time] .% 1.0 .== 0, :] # Not(:time)]

    ## Compute norm of difference
    # if plot_comparison # Code snippte to visualize differences for development purposes
        # using Plots
        # pl = plot(sim_θ[Not(1),:], line = :solid,
        #     labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
        # plot!(ref_θ, line = :dash, color = :black,
        #     labels = "LWFBrook90R: " .* string.(depth_to_read_out_mm) .* " mm")
        # plot(HydrusSolution_sparseTime_θ_wide[:,:time],
        #     Matrix(HydrusSolution_sparseTime_θ_wide[:,Not(:time)]),
        #     line = :dash, color = :green)
        # plot!(HydrusSolution_denseTime_θ_wide[:,:time],
        #     Matrix(HydrusSolution_denseTime_θ_wide[:,Not(:time)]),
        #     line = :dash, color = :red)
        # plot(HydrusSolution_sparseTime_δ2H_wide[:,:time],
        #     Matrix(HydrusSolution_sparseTime_δ2H_wide[:,Not(:time)]),
        #     line = :dash, color = :red)
        # plot(HydrusSolution_sparseTime_δ18O_wide[:,:time],
        #     Matrix(HydrusSolution_sparseTime_δ18O_wide[:,Not(:time)]),
        #     line = :dash, color = :red)
    # end

    # TODO(bernhard): Check why currently sim_θ needs to be shifted by 1 day.
    #                 Is it the initial conditions? (That are not reported by LWFBrook90R?)
    return (
        sim = (θ = sim_θ[Not(1),:],   ψ = sim_ψ[Not(1),:],   δ18O = sim_δ18O[Not(1),:],   δ2H = sim_δ2H[Not(1),:]),
        # also remove from ref_θ a day at the end to have same number:
        ref = (θ = ref_θ[Not(end),:], ψ = ref_ψ[Not(end),:], δ18O = ref_δ18O[Not(end),:], δ2H = ref_δ2H[Not(end),:]),
        hyd = (θ = hyd_θ,             ψ = hyd_ψ,             δ18O = hyd_δ18O,             δ2H = hyd_δ2H,
               θdense = HydrusSolution_denseTime_θ_wide) # this is densely output in time
    )
end

function prepare_sim_and_ref_for_BEA_2016(
    folder_with_sim_input_and_ref_output,
    input_prefix)

    # folder_with_sim_input_and_ref_output = "test-assets/BEA-2016"
    # input_prefix = "BEA2016-reset-FALSE"

    # Run  simulation
    sim_sol, _, _ = run_simulation(
    [joinpath(folder_with_sim_input_and_ref_output, "input-files/"),
    input_prefix]
    )

    # Postprocess simulation
    ## Aboveground

    ## Belowground
    ### Extract certain depths of simulation
    depth_to_read_out_mm = [100 500 1000 1500 1900]
    idx_sim_sol = find_indices(depth_to_read_out_mm, sim_sol)

    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
            get_auxiliary_variables(sim_sol)

    sim = (θ = insertcols!(DataFrame(u_aux_θ[:,idx_sim_sol], :auto),
                :time => sim_sol.t),
           ψ = insertcols!(DataFrame(u_aux_PSIM[:,idx_sim_sol], :auto),
                :time => sim_sol.t),
           above = insertcols!(DataFrame([sim_sol.u[t_idx][u_idx] for
                    t_idx in 1:length(sim_sol),
                    u_idx in [1,2,3,4,5,6]],
                ["GWAT","INTS","INTR", "SNOW", "CC", "SNOWLQ"]),
            :time => sim_sol.t))

    # Load LWFbrook90R solution
    ## Aboveground
    ref_above_1 = read_LWFBrook90R_dailyCSV(path =
        joinpath(
            folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
            input_prefix*"_NLAYER7_LWFBrook90R-0.4.5daily_output.csv"))
    ref_above_2 = read_LWFBrook90R_dailyCSV(path =
        joinpath(
            folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
            input_prefix*"_NLAYER14_LWFBrook90R-0.4.5daily_output.csv"))
    ref_above_3 = read_LWFBrook90R_dailyCSV(path =
        joinpath(
            folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
            input_prefix*"_NLAYER21_LWFBrook90R-0.4.5daily_output.csv"))
    ref_above_4 = read_LWFBrook90R_dailyCSV(path =
        joinpath(
            folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
            input_prefix*"_NLAYER70_LWFBrook90R-0.4.5daily_output.csv"))

    ## Belowground
    ref_below_1 = read_LWFBrook90R_layerCSV_extract_depths(;
        path = joinpath(folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
                        input_prefix*"_NLAYER7_LWFBrook90R-0.4.5layer_output.csv"),
        depth_to_read_out_mm = [100 500 1000 1500 1900])

    ref_below_2 = read_LWFBrook90R_layerCSV_extract_depths(;
        path = joinpath(folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
                        input_prefix*"_NLAYER14_LWFBrook90R-0.4.5layer_output.csv"),
        depth_to_read_out_mm = [100 500 1000 1500 1900])

    ref_below_3 = read_LWFBrook90R_layerCSV_extract_depths(;
        path = joinpath(folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
                        input_prefix*"_NLAYER21_LWFBrook90R-0.4.5layer_output.csv"),
        depth_to_read_out_mm = [100 500 1000 1500 1900])

    ref_below_4 = read_LWFBrook90R_layerCSV_extract_depths(;
        path = joinpath(folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
                        input_prefix*"_NLAYER70_LWFBrook90R-0.4.5layer_output.csv"),
        depth_to_read_out_mm = [100 500 1000 1500 1900])

    # Make simulation and reference comparable (formatting)
    # TODO(bernhard): currently sim_θ needs to be shifted by 1 day.
    #                 is it the initial conditions? (That are not reported by LWFBrook90R?)
    return (
        (θ = sim.θ[Not(1),:], ψ = sim.ψ[Not(1),:], above = sim.above),
        # also remove from ref_below_.XXX a day at the end to have same number
        (θ = ref_below_1.θ[Not(end),:], ψ = ref_below_1.ψ[Not(end),:], above = ref_above_1),
        (θ = ref_below_2.θ[Not(end),:], ψ = ref_below_2.ψ[Not(end),:], above = ref_above_2),
        (θ = ref_below_3.θ[Not(end),:], ψ = ref_below_3.ψ[Not(end),:], above = ref_above_3),
        (θ = ref_below_4.θ[Not(end),:], ψ = ref_below_4.ψ[Not(end),:], above = ref_above_4)
    )

    # # Code snippte to visualize differences for development purposes
    # using Plots
    # ## Aboveground
    # pl_ab_all = plot(sim_sol; vars = [1, 2, 3, 4, 5, 6],
    #      label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])
    # # alternatively:
    # plot(sim.above.time,
    #      Matrix(sim.above[:,Not(:time)]),
    #      label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])

    # pl_ab_3 = plot(sim_sol; vars = [2, 3, 4], label=["INTS (mm)" "INTR (mm)" "SNOW (mm)"])
    # pl_ab_2 = plot(sim_sol; vars = [2, 3],    label=["INTS (mm)" "INTR (mm)"])

    # # add LWFBrook90References
    # plot!(pl_ab_all,
    #     [ref_above_1.gwat,
    #       ref_above_1.intr,
    #       ref_above_1.ints,
    #       ref_above_1.snow], label = "LWFBrook90R", color = :black, line = :dash)
    # plot!(pl_ab_3,
    #     [ref_above_1.intr,
    #       ref_above_1.ints,
    #       ref_above_1.snow], label = "LWFBrook90R", color = :black, line = :dash)
    # plot!(pl_ab_2,
    #     [ref_above_1.intr,
    #       ref_above_1.ints], label = "LWFBrook90R", color = :black, line = :dash)

    # ## Belowground
    # ## θ
    # pl = plot(Matrix(sim.θ[Not(1),Not(:time)]), line = :solid,
    #     labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    # plot!(ref_below_1.θ, line = :dash, color = :black,
    #     labels = "LWFBrook90R: " .* string.(depth_to_read_out_mm) .* " mm")
    # plot!(ref_below_2.θ, line = :dash, color = :black,
    #     labels = "LWFBrook90R: " .* string.(depth_to_read_out_mm) .* " mm")
    # plot!(ref_below_3.θ, line = :dash, color = :black,
    #     labels = "LWFBrook90R: " .* string.(depth_to_read_out_mm) .* " mm")
    # plot!(ref_below_4.θ, line = :dash, color = :black,
    #     labels = "LWFBrook90R: " .* string.(depth_to_read_out_mm) .* " mm")
    # ## ψ
    # pl = plot(Matrix(sim.ψ[Not(1),Not(:time)]), line = :solid,
    #     labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    # plot!(ref_below_1.ψ, line = :dash, color = :black,
    #     labels = "LWFBrook90R: " .* string.(depth_to_read_out_mm) .* " mm")
    # plot!(ref_below_2.ψ, line = :dash, color = :black,
    #     labels = "LWFBrook90R: " .* string.(depth_to_read_out_mm) .* " mm")
    # plot!(ref_below_3.ψ, line = :dash, color = :black,
    #     labels = "LWFBrook90R: " .* string.(depth_to_read_out_mm) .* " mm")
    # plot!(ref_below_4.ψ, line = :dash, color = :black,
    #     labels = "LWFBrook90R: " .* string.(depth_to_read_out_mm) .* " mm")

end


function read_LWFBrook90R_dailyCSV(;path)
    # path = joinpath(
    #     folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
    #     input_prefix*"_NLAYER7_LWFBrook90R-0.4.5daily_output.csv")

    ## Read in csv
    referenceSolution_daily = DataFrame(File(path))

    ref_aboveground = referenceSolution_daily[:,
            [#:yr, :mo, :da,
             :doy, :intr, :ints, :snow, :gwat,
                #:snowlq, :cc
            ]]
    rename!(ref_aboveground, :doy => :time)
end

function read_LWFBrook90R_layerCSV_extract_depths(;path, depth_to_read_out_mm)
    ## Read in csv
    referenceSolution_layer = DataFrame(File(path))

    ## Extract certain depths
    idx_referenceSolution = []
    for curr_depth_mm in
        # Only read out values that are within the simulation domain
        depth_to_read_out_mm[depth_to_read_out_mm .<= maximum(-referenceSolution_layer.lower)]

        # idx_to_read_out = findfirst(curr_depth_mm .<= y)
        append!(idx_referenceSolution, findfirst(curr_depth_mm .<= -referenceSolution_layer.lower))
    end

    ## Process into wide data frame
    ref_θ = unstack(
        referenceSolution_layer[:,[:yr, :mo, :da, :doy, :nl, :theta]],
        [:yr, :mo, :da, :doy], # ID variable
        :nl,   # variable that will be spread out into column names
        :theta,# value that will be filled into the wide table
        renamecols=x->Symbol(:nl_, x)
        )[:,
            [:doy; Symbol.("nl_" .* string.(idx_referenceSolution))]]
            # Symbol.("nl_" .* string.(idx_referenceSolution))]
            # [:yr; :mo; :da; :doy; Symbol.("nl_" .* string.(idx_referenceSolution))]])
    rename!(ref_θ, :doy => :time)

    ref_ψ = unstack(
            referenceSolution_layer[:,[:yr, :mo, :da, :doy, :nl, :psimi]],
            [:yr, :mo, :da, :doy], # ID variable
            :nl,   # variable that will be spread out into column names
            :psimi,# value that will be filled into the wide table
            renamecols=x->Symbol(:nl_, x)
            )[:,
                [:doy; Symbol.("nl_" .* string.(idx_referenceSolution))]]
                # Symbol.("nl_" .* string.(idx_referenceSolution))]
                # [:yr; :mo; :da; :doy; Symbol.("nl_" .* string.(idx_referenceSolution))]]
    rename!(ref_ψ, :doy => :time)

    return (θ = ref_θ, ψ = ref_ψ)
end
