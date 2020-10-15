using Pkg.TOML, StaticArrays

@testset "MagneticFields.jl" begin
    t = randn()
    r = rand(SVector{3,Float64})
    G = SMatrix{3,3,Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -2.0])

    params = TOML.parsefile("MagneticFields.toml")
	
	staticHomogeneousField = MPISimulations.initialize(MPISimulations.StaticHomogeneousField,"StaticHomogeneousField",params,Float64)
    @test MPISimulations.magneticFieldStrength(staticHomogeneousField,r) ≈ SVector{3,Float64}(0.12,0.12,0.12)
    @test MPISimulations.magneticFieldStrength(staticHomogeneousField,r,t) ≈ SVector{3,Float64}(0.12,0.12,0.12)

    staticGradientField = MPISimulations.initialize(MPISimulations.StaticGradientField,"StaticGradientField",params,Float64)
    @test MPISimulations.magneticFieldStrength(staticGradientField,r) ≈ G*r
    @test MPISimulations.magneticFieldStrength(staticGradientField,r,t) ≈ G*r

    sinusoidalField1 = MPISimulations.initialize(MPISimulations.SinusoidalField,"SinusoidalField",params,Float64)
    @test MPISimulations.magneticFieldStrength(sinusoidalField1,r,t) ≈ sin(2*pi*25e3*t+π)*SVector{3,Float64}(0.12,0.12,0.12)

    spatialFieldProfile = MPISimulations.initialize(MPISimulations.StaticHomogeneousField,"SinusoidalField.spatialFieldProfile.StaticHomogeneousField",params,Float64)
    sinusoidalField2 = MPISimulations.SinusoidalField(spatialFieldProfile,2.5e6,100,1.0)
    @test sinusoidalField1.baseFrequency == sinusoidalField2.baseFrequency
    @test sinusoidalField1.divider == sinusoidalField2.divider
    @test sinusoidalField1.doubleFrequency == sinusoidalField2.doubleFrequency
    @test sinusoidalField1.relativePhase == sinusoidalField2.relativePhase
    @test sinusoidalField1.spatialFieldProfile.magneticFieldStrength == sinusoidalField2.spatialFieldProfile.magneticFieldStrength
end