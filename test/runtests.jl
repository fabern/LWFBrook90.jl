# fabian.bernhard@wsl.ch, 2021-02-01

using LWFBrook90Julia
using Test

@testset "LWFBrook90Julia.jl" begin
    # Write your tests here.

    # 2x + 3y
    @test my_f(2,1) == 7
    @test my_f(2,3) == 13
end
