# Integration tests

# - _Unit testing_ asserts that individual pieces of a project work as expected. (developers
#       perspective)
# - _Integration testing_ asserts that they fit together as expected. Also known as
#       _functional tests_, they cover entire use cases (user perspective). For LWFBrook90.jl
#       these are tests that are compared to e.g. LWFBrook90R or Hydrus.
# - _Regression testing_ asserts that behavior is unchanged over time. Also known as
#       _reference tests_.

include("fct-helpers-for-integration-tests.jl")

# TODO(bernhard): include further quantities such as ψ etc...
# TODO(bernhard): include further "true" values such as the Hammel-2001 fig 2.4/2.5 infiltration test

# NOTE: locally, i.e. not on CI system, one might need to do manually cd("test")
@testset "Hammel-2001-θ" begin

    @show pwd() # This is to help get the folder right.

    # Check the RMSE of θ in simulations is below a limit
    @githash_time sim1, ref1, hyd1 = prepare_θ_from_sim_and_reference("test-assets/Hammel-2001","Hammel_loam-NLayer-27-RESET=TRUE")
        # amberMBP-git-c4275ee: 4.394856 seconds (19.52 M allocations: 2.898 GiB, 16.30% gc time)
    @githash_time sim2, ref2, hyd2 = prepare_θ_from_sim_and_reference("test-assets/Hammel-2001","Hammel_loam-NLayer-103-RESET=TRUE")
        # amberMBP-git-c4275ee: 13.704865 seconds (19.92 M allocations: 9.220 GiB, 15.41% gc time, 0.23% compilation time)
    # @githash_time sim3, ref3, hyd3 = prepare_θ_from_sim_and_reference("test-assets/Hammel-2001","Hammel_loam-NLayer-400-RESET=TRUE")
        # not run

    @githash_time sim4, ref4, hyd4 = prepare_θ_from_sim_and_reference("test-assets/Hammel-2001","Hammel_sand-NLayer-27-RESET=TRUE")
        # amberMBP-git-c4275ee: 0.913546 seconds (4.04 M allocations: 605.712 MiB, 14.49% gc time, 3.67% compilation time)
    @githash_time sim5, ref5, hyd5 = prepare_θ_from_sim_and_reference("test-assets/Hammel-2001","Hammel_sand-NLayer-103-RESET=TRUE")
        # amberMBP-git-c4275ee: 2.788329 seconds (3.92 M allocations: 1.817 GiB, 15.26% gc time)
    # @githash_time sim6, ref6, hyd6 = prepare_θ_from_sim_and_reference("test-assets/Hammel-2001","Hammel_sand-NLayer-400-RESET=TRUE")
        # not run

    function RMS_differences(sim, ref)
        # computes root mean squared differences
        differences = sim .- ref
        mean_squared = sum(differences.^2) / length(differences)

        return sqrt(mean_squared)
    end

    # Use sensible accuracy values to compare the two solutions (e.g. θ of 0.02, and ψ of 1 kPa)
    # (e.g. RMSE(reference, simulated) < [hardcoded_value])
    # Floating point accuracy is used by regression tests.
    # (Or could be used, once I compare with a very high-fidelity solution, of which I known
    #  that my code should converge to. E.g. a high-resolution Hydrus simulation.)

    # Compare with LWFBrook90R as reference solution
    @test RMS_differences(sim1, ref1) < 0.03
    @test RMS_differences(sim2, ref2) < 0.015
    # @test RMS_differences(sim3, ref3) < 0.0015
    @test RMS_differences(sim4, ref4) < 0.02
    @test RMS_differences(sim5, ref5) < 0.01
    # @test RMS_differences(sim6, ref6) < 0.005

    # Compare with Hydrus1D
    @test RMS_differences(sim1[Not(end),:], hyd1) < 0.005
    @test RMS_differences(sim2[Not(end),:], hyd2) < 0.002
    # @test RMS_differences(sim3[Not(end),:], hyd3) < 0.001
    @test RMS_differences(sim4[Not(end),:], hyd4) < 0.006
    @test RMS_differences(sim5[Not(end),:], hyd5) < 0.007
    # @test RMS_differences(sim6[Not(end),:], hyd6) < 0.005

    # if some error appears, the following code can be used to plot the solutions
    # using Plots
    # depth_to_read_out_mm = [100 500 1000 1500 1900]
    # pl = plot(sim4, line = :solid, labels = "LWFBrook90.jl: " .* string.(depth_to_read_out_mm) .* " mm")
    # plot!(ref4, line = :dash, color = :black, labels = "LWFBrook90R: " .* string.(depth_to_read_out_mm) .* " mm")
    # plot!(hyd4, line = :dash, color = :green, labels = "Hydrus: " .* string.(depth_to_read_out_mm) .* " mm")

end
