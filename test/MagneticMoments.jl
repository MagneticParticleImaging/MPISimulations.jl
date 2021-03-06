using Pkg.TOML, StaticArrays, LinearAlgebra

@testset "MagneticMoments.jl" begin
	params = TOML.parsefile("MagneticMoments.toml")
	
	magneticMoment = MPISimulations.initialize(MPISimulations.Langevin,"Langevin",params,Float64)
	for T in [Float64,Float32,]
		magneticMomentT = MPISimulations.Langevin(elementType=T)
		@test magneticMomentT.msat ≈ magneticMoment.msat
		@test magneticMomentT.beta ≈ magneticMoment.beta
	end
	
	t = randn()
	r = rand(SVector{3,Float64})
	H = rand(SVector{3,Float64})
	H = zero(SVector{3,Float64})
	@test MPISimulations.meanMagneticMoment(magneticMoment,H,r,t) ≈ zero(SVector{3,Float64})
	H = SVector{3,Float64}(1e5,0,0)
	@test MPISimulations.meanMagneticMoment(magneticMoment,H,r,t) ≈ SVector{3,Float64}(magneticMoment.msat,0,0)
	H = SVector{3,Float64}(0,1e5,0)
	@test MPISimulations.meanMagneticMoment(magneticMoment,H,r,t) ≈ SVector{3,Float64}(0,magneticMoment.msat,0)
	H = SVector{3,Float64}(0,0,1e5)
	@test MPISimulations.meanMagneticMoment(magneticMoment,H,r,t) ≈ SVector{3,Float64}(0,0,magneticMoment.msat)
	H = randn(SVector{3,Float64})
	x = magneticMoment.beta*norm(H)
	@test MPISimulations.meanMagneticMoment(magneticMoment,H,r,t) ≈ magneticMoment.msat*(coth(x) - 1/x)*normalize(H)
end
