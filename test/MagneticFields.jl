using Pkg.TOML, StaticArrays

@testset "MagneticFields.jl" begin
    t = randn()
    r = rand(SVector{3,Float64})
    G = SMatrix{3,3,Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -2.0])

    params = TOML.parsefile("MagneticFields.toml")
	
	staticHomogeneousField = MPISimulations.initialize(MPISimulations.StaticHomogeneousField,params,Float64)
    @test MPISimulations.magneticFieldStrength(staticHomogeneousField,r) ≈ SVector{3,Float64}(0.12,0.12,0.12)
    @test MPISimulations.magneticFieldStrength(staticHomogeneousField,r,t) ≈ SVector{3,Float64}(0.12,0.12,0.12)

    staticGradientField = MPISimulations.initialize(MPISimulations.StaticGradientField,params,Float64)
    @test MPISimulations.magneticFieldStrength(staticGradientField,r) ≈ G*r
    @test MPISimulations.magneticFieldStrength(staticGradientField,r,t) ≈ G*r

    sinusoidalExcitationField = MPISimulations.SinusoidalExcitationField(staticHomogeneousField,25e3,0.0)
    @test MPISimulations.magneticFieldStrength(sinusoidalExcitationField,r,t) ≈ sin(2*pi*25e3*t)*SVector{3,Float64}(0.12,0.12,0.12)
end