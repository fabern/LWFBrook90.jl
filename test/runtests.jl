# fabian.bernhard@wsl.ch, 2021-02-01

using LWFBrook90Julia
using Test

# @testset "Module WAT" begin
@testset "KKMEAN" begin
    @test LWFBrook90Julia.WAT.KKMEAN(50.,  50., 1., 1.) ≈ 50
    @test LWFBrook90Julia.WAT.KKMEAN(50.,  50., 1., 5.) ≈ 50
    @test LWFBrook90Julia.WAT.KKMEAN(10., 100., 1., 1.) ≈ 31.6227766
    @test LWFBrook90Julia.WAT.KKMEAN(10., 100., 1., 5.) ≈ 14.6779926
    @test LWFBrook90Julia.WAT.KKMEAN(10., 100., 5., 1.) ≈ 68.1292069
end