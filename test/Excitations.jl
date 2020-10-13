using Pkg.TOML, StaticArrays

@testset "Excitations.jl" begin
    params = TOML.parsefile("Excitations.toml")
	
    excitationField1D = MPISimulations.initialize(MPISimulations.ExcitationField1D,"ExcitationField1D",params,Float64)
    r = zero(SVector{3,Float64})
    @test MPISimulations.magneticExcitationFieldSequence(excitationField1D,r) == SVector{3,Float64}[SVector{3,Float64}(0.12, 0.0, 0.0),SVector{3,Float64}(0.0, 0.0, 0.0),SVector{3,Float64}(-0.12, 0.0, 0.0)]
    r = SVector{3,Float64}(0.0,0.0,0.1)
    @test MPISimulations.magneticExcitationFieldSequence(excitationField1D,r) == SVector{3,Float64}[SVector{3,Float64}(0.12, 0.0, -0.2),SVector{3,Float64}(0.0, 0.0, -0.2),SVector{3,Float64}(-0.12, 0.0, -0.2)]

    staticGradientField = MPISimulations.initialize(MPISimulations.StaticGradientField,"StaticGradientField",params,Float64)
    sinusoidalField = MPISimulations.initialize(MPISimulations.SinusoidalField,"SinusoidalField",params,Float64)
    timepoints = range(0, 1/25e3, length=5)
    excitationField1D = MPISimulations.ExcitationField1D(staticGradientField,sinusoidalField,timepoints)
    @test excitationField1D.selectionField == staticGradientField
    @test excitationField1D.excitationField == sinusoidalField.spatialFieldAmplitude
    @test excitationField1D.excitationStrengths ≈ Float64[0,1,0,-1,0]
end