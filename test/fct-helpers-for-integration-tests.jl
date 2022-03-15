############################################################################################
############################################################################################
############################################################################################
# Helpers for integration testing

function prepare_θ_from_sim_and_reference(
    folder_with_sim_input_and_ref_output,
    input_prefix)

    #e.g. folder_with_sim_input_and_ref_output = "test/test-assets/Hammel-2001"
    #e.g. input_prefix = "Hammel_loam-NLayer-27-RESET=TRUE"

    # Run  simulation
    sim_sol, _, _ = run_simulation(
    [joinpath(folder_with_sim_input_and_ref_output, "input-files/"),
    input_prefix] # "Hammel_loam-NLayer-27-RESET=TRUE"]
    )

    # Load reference solution
    ## Load Hydrus solution
    HydrusSolution = DataFrame(
        File(
            joinpath(
                "test-assets/Hammel-2001", "output_Hydrus1D",
            "Hammel_Test_"*uppercasefirst(replace(input_prefix, r"Hammel_([[:alnum:]]*)-NLayer.*" => s"\g<1>")),
            "", "Obs_Node_processed.csv")
            ))
    ## Load LWFbrook90R solution
    # referenceSolution_daily = DataFrame(
    #     File(joinpath(
    #         folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
    #         input_prefix*"_OUTPUT-LWFBrook90R-0.4.5-daily_output.csv")))
    referenceSolution_layer = DataFrame(
        File(joinpath(
            folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
            input_prefix*"_OUTPUT-LWFBrook90R-0.4.5-layer_output.csv")))
    # SOME NOTES:
    # z_lower_mm = round.(-cumsum(sim_sol.prob.p[1][1].p_THICK); digits=5)
    # z_upper_mm = [0; z_lower_mm[Not(end)]]
    # z_midpoints_mm = (z_upper_mm .+ z_lower_mm)./2

    # Compare simulation and reference
    ## Extract certain depths
    depth_to_read_out_mm = [100 500 1000 1500 1900]
    idx = find_indices(depth_to_read_out_mm, sim_sol)
    ### Ref
    ref_θ = Matrix(unstack(
        referenceSolution_layer[:,[:yr, :mo, :da, :doy, :nl, :theta]],
        [:yr, :mo, :da, :doy], # ID variable
        :nl,   # variable that will be spread out into column names
        :theta,# value that will be filled into the wide table
        renamecols=x->Symbol(:nl_, x)
        )[:,idx])
    ### Sim
    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
            get_auxiliary_variables(sim_sol)
    sim_θ = u_aux_θ[:,idx]

    ### Hydrus
    HydrusSolution_wide = unstack(
        HydrusSolution[:, [:time, :depth, :theta]],
        [:time], # ID variable
        :depth,  # variable that will be spread out into column names
        :theta,  # value that will be filled into the wide table
        renamecols=x->Symbol(:depth_, x)
        )

    HydrusSolution_final =
        Matrix(HydrusSolution_wide[HydrusSolution_wide[:, :time] .% 1.0 .== 0, Not(:time)])

    ## Compute norm of difference
    # if plot_comparison # Code snippte to visualize differences for development purposes
    #     pl = plot(sim_θ[Not(1),:], line = :solid,
    #         labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    #     plot!(ref_θ, line = :dash, color = :black,
    #         labels = "LWFBrook90R: " .* string.(depth_to_read_out_mm) .* " mm")
    # end

    # TODO(bernhard): currently sim_θ needs to be shifted by 1 day.
    #                 is it the initial conditions? (That are not reported by LWFBrook90R?)
    return (
        sim_θ[Not(1),:],
        ref_θ[Not(end),:], # also remove from ref_θ a day at the end to have same number
        HydrusSolution_final
    )
end

function prepare_sim_and_ref_for_BEA_2016(
    folder_with_sim_input_and_ref_output,
    input_prefix)

    #e.g. folder_with_sim_input_and_ref_output = "test/test-assets/Hammel-2001"
    #e.g. input_prefix = "Hammel_loam-NLayer-27-RESET=TRUE"

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

    sim = (θ = u_aux_θ[:,idx_sim_sol],
           ψ = u_aux_PSIM[:,idx_sim_sol],
           above = [sim_sol.u[t_idx][u_idx] for
                    t_idx in 1:length(sim_sol),
                    u_idx in [1,2,3,4,5,6]])

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
    ref_below_1 = read_LWFBrook90R_layerCSV_extract_theta_depths(;
        path = joinpath(folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
                        input_prefix*"_NLAYER7_LWFBrook90R-0.4.5layer_output.csv"),
        depth_to_read_out_mm = [100 500 1000 1500 1900])

    ref_below_2 = read_LWFBrook90R_layerCSV_extract_theta_depths(;
        path = joinpath(folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
                        input_prefix*"_NLAYER14_LWFBrook90R-0.4.5layer_output.csv"),
        depth_to_read_out_mm = [100 500 1000 1500 1900])

    ref_below_3 = read_LWFBrook90R_layerCSV_extract_theta_depths(;
        path = joinpath(folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
                        input_prefix*"_NLAYER21_LWFBrook90R-0.4.5layer_output.csv"),
        depth_to_read_out_mm = [100 500 1000 1500 1900])

    ref_below_4 = read_LWFBrook90R_layerCSV_extract_theta_depths(;
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
    # ## Aboveground
    # pl_ab_all = plot(sim_sol; vars = [1, 2, 3, 4, 5, 6],
    #      label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])
    # # alternatively:
    # plot(sim.above,
    #     label=["GWAT (mm)" "INTS (mm)" "INTR (mm)" "SNOW (mm)" "CC (MJ/m2)" "SNOWLQ (mm)"])
    #
    # pl_ab_3 = plot(sim_sol; vars = [2, 3, 4], label=["INTS (mm)" "INTR (mm)" "SNOW (mm)"])
    # pl_ab_2 = plot(sim_sol; vars = [2, 3],    label=["INTS (mm)" "INTR (mm)"])
    #
    # # add LWFBrook90References
    # plot!(pl_ab_all,
    #     [ref_above_1.gwat,
    #       ref_above_1.intr,
    #       ref_above_1.ints,
    #       ref_above_1.snow], label = "LWFBrook90R", line = :dash)
    # plot!(pl_ab_3,
    #     [ref_above_1.intr,
    #       ref_above_1.ints,
    #       ref_above_1.snow], label = "LWFBrook90R", line = :dash)
    # plot!(pl_ab_2,
    #     [ref_above_1.intr,
    #       ref_above_1.ints], label = "LWFBrook90R", line = :dash)
    #
    # ## Belowground
    # ## θ
    # pl = plot(sim.θ[Not(1),:], line = :solid,
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
    # pl = plot(sim.ψ[Not(1),:], line = :solid,
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
            [:yr, :mo, :da, :doy, :intr, :ints, :snow, :gwat,
                #:snowlq, :cc
            ]]
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
    ref_θ = Matrix(unstack(
        referenceSolution_layer[:,[:yr, :mo, :da, :doy, :nl, :theta]],
        [:yr, :mo, :da, :doy], # ID variable
        :nl,   # variable that will be spread out into column names
        :theta,# value that will be filled into the wide table
        renamecols=x->Symbol(:nl_, x)
        )[:,
            Symbol.("nl_" .* string.(idx_referenceSolution))])
            # [:yr; :mo; :da; :doy; Symbol.("nl_" .* string.(idx_referenceSolution))]]))

    ref_ψ = Matrix(unstack(
            referenceSolution_layer[:,[:yr, :mo, :da, :doy, :nl, :psimi]],
            [:yr, :mo, :da, :doy], # ID variable
            :nl,   # variable that will be spread out into column names
            :psimi,# value that will be filled into the wide table
            renamecols=x->Symbol(:nl_, x)
            )[:,
                Symbol.("nl_" .* string.(idx_referenceSolution))])
                # [:yr; :mo; :da; :doy; Symbol.("nl_" .* string.(idx_referenceSolution))]]))

    return (θ = ref_θ, ψ = ref_ψ)
end

