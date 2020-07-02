using Pkg.TOML, StaticArrays, LinearAlgebra

@testset "TracerDistributions.jl" begin
	params = TOML.parsefile("TracerDistributions.toml")

	aob = MPISimulations.initialize(MPISimulations.AxisOrientedBox,params,Float64)
	sphere = MPISimulations.initialize(MPISimulations.Sphere,params,Float64)
	t = randn()
	r = zeros(SVector{3,Float64})
	@test MPISimulations.isInside(aob,r)
	@test MPISimulations.isInside(sphere,r)
	r = SVector{3,Float64}(1e-3,0,0)
	@test MPISimulations.isInside(aob,r)
	@test MPISimulations.isInside(sphere,r)
	r = SVector{3,Float64}(1.01e-3,0,0)
	@test !MPISimulations.isInside(aob,r)
	@test !MPISimulations.isInside(sphere,r)

	r = zeros(SVector{3,Float64})
	@test MPISimulations.tracerConcentration(aob,r,t) == 0.42
	@test MPISimulations.tracerConcentration(sphere,r,t) == 0.42
	r = SVector{3,Float64}(1.01e-3,0,0)
	@test MPISimulations.tracerConcentration(aob,r,t) == 0.0
	@test MPISimulations.tracerConcentration(sphere,r,t) == 0.0
end
