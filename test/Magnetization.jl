using Pkg.TOML, StaticArrays, LinearAlgebra

@testset "Magnetization.jl" begin
    params = TOML.parsefile("Magnetization.toml")
	
    separableMagnetization = MPISimulations.initialize(MPISimulations.SeparableMagnetization,"SeparableMagnetization",params,Float64)
    r = zeros(SVector{3,Float64})
    t = randn()
    H = SVector{3,Float64}(1,0,0)
    x = 2e3*norm(H)
    
    @test !MPISimulations.isZero(separableMagnetization,r)
    @test MPISimulations.magnetizationField(separableMagnetization,H,r,t) ≈ 42e-2*7e-18*(coth(x) - 1/x)*normalize(H)
	r = SVector{3,Float64}(1e-3,0,0)
    @test !MPISimulations.isZero(separableMagnetization,r)
    @test MPISimulations.magnetizationField(separableMagnetization,H,r,t) ≈ 42e-2*7e-18*(coth(x) - 1/x)*normalize(H)
	r = SVector{3,Float64}(1.01e-3,0,0)
    @test MPISimulations.isZero(separableMagnetization,r)
    @test MPISimulations.magnetizationField(separableMagnetization,H,r,t) ≈ zeros(SVector{3,Float64})
end