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
    # reference_loam_27_daily = DataFrame(
    #     File(joinpath(
    #         folder_with_sim_input_and_ref_output, "output_LWFBrook90R",
    #         input_prefix*"_OUTPUT-LWFBrook90R-0.4.5-daily_output.csv")))
    reference_loam_27_layer = DataFrame(
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
        reference_loam_27_layer[:,[:yr, :mo, :da, :doy, :nl, :theta]],
        [:yr, :mo, :da, :doy], # ID variable
        :nl,   # variable that will be spread out into column names
        :theta,# value that will be filled into the wide table
        renamecols=x->Symbol(:nl_, x)
        )[:,idx])
    ### Sim
    (u_SWATI, u_aux_WETNES, u_aux_PSIM, u_aux_PSITI, u_aux_θ, p_fu_KK) =
            get_auxiliary_variables(sim_sol)
    sim_θ = u_aux_θ[:,idx]

    ## Compute norm of difference
    # if plot_comparison # Code snippte to visualize differences for development purposes
    #     pl = plot(sim_θ[Not(1),:], line = :solid,
    #         labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    #     plot!(ref_θ, line = :dash, color = :black,
    #         labels = "LWFBrook90R: " .* string.(depth_to_read_out_mm) .* " mm")
    # end

    # TODO(bernhard): currently sim_θ needs to be shifted by 1 day.
    #                 is it the initial conditions? (That are not reported by LWFBrook90R?)
    return (sim_θ[Not(1),:], ref_θ[Not(end),:]) # also remove from ref_θ a day at the end to have same number
end