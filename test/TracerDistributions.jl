using Pkg.TOML, StaticArrays

@testset "TracerDistributions.jl" begin
	params = TOML.parsefile("TracerDistributions.toml")

	aob = MPISimulations.initialize(MPISimulations.AxisOrientedBox,"AxisOrientedBox",params,Float64)
	sphere = MPISimulations.initialize(MPISimulations.Sphere,"Sphere",params,Float64)
	t = randn()
	r = zeros(SVector{3,Float64})
	@test !MPISimulations.isZero(aob,r)
	@test !MPISimulations.isZero(sphere,r)
	r = SVector{3,Float64}(1e-3,0,0)
	@test !MPISimulations.isZero(aob,r)
	@test !MPISimulations.isZero(sphere,r)
	r = SVector{3,Float64}(1.01e-3,0,0)
	@test MPISimulations.isZero(aob,r)
	@test MPISimulations.isZero(sphere,r)

	r = zeros(SVector{3,Float64})
	@test MPISimulations.tracerConcentration(aob,r,t) == 1.8141555811322022e15
	@test MPISimulations.tracerConcentration(sphere,r,t) == 1.8141555811322022e15
	r = SVector{3,Float64}(1.01e-3,0,0)
	@test MPISimulations.tracerConcentration(aob,r,t) == 0.0
	@test MPISimulations.tracerConcentration(sphere,r,t) == 0.0

	#Misc
	@test MPISimulations.numberOfParticles(1e-3,25e-9) == 1.8141352092994865e15
end
