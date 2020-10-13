using Pkg.TOML, StaticArrays

@testset "Excitations.jl" begin
    params = TOML.parsefile("Excitations.toml")
	
    excitationField1D = MPISimulations.initialize(MPISimulations.ExcitationField1D,"ExcitationField1D",params,Float64)
    r = zero(SVector{3,Float64})
    @test MPISimulations.magneticExcitationFieldSequence(excitationField1D,r) == SVector{3,Float64}[SVector{3,Float64}(0.12, 0.0, 0.0),SVector{3,Float64}(0.0, 0.0, 0.0),SVector{3,Float64}(-0.12, 0.0, 0.0)]
    r = SVector{3,Float64}(0.0,0.0,0.1)
    @test MPISimulations.magneticExcitationFieldSequence(excitationField1D,r) == SVector{3,Float64}[SVector{3,Float64}(0.12, 0.0, -0.2),SVector{3,Float64}(0.0, 0.0, -0.2),SVector{3,Float64}(-0.12, 0.0, -0.2)]

    # create excitation from static and sinusoidal field
    staticGradientField = MPISimulations.initialize(MPISimulations.StaticGradientField,"StaticGradientField",params,Float64)
    sinusoidalField1 = MPISimulations.initialize(MPISimulations.SinusoidalField,"SinusoidalField1",params,Float64)
    sinusoidalField2 = MPISimulations.initialize(MPISimulations.SinusoidalField,"SinusoidalField2",params,Float64)
    sinusoidalField3 = MPISimulations.initialize(MPISimulations.SinusoidalField,"SinusoidalField3",params,Float64)
    timepoints = range(0, 1, length=10)

    excitationField1D = MPISimulations.ExcitationField(staticGradientField,sinusoidalField1,timepoints)
    @test excitationField1D.selectionField == staticGradientField
    @test excitationField1D.excitationField == sinusoidalField1.spatialFieldAmplitude
    @test excitationField1D.excitationStrengths ≈ map(t->sin(2*pi*24509.803921568626*t+1),timepoints)

    excitationField2D = MPISimulations.ExcitationField(staticGradientField,sinusoidalField1,sinusoidalField2,timepoints)
    @test excitationField2D.selectionField == staticGradientField
    @test excitationField2D.excitationField1 == sinusoidalField1.spatialFieldAmplitude
    @test excitationField2D.excitationField2 == sinusoidalField2.spatialFieldAmplitude
    @test excitationField2D.excitationStrengths1 ≈ map(t->sin(2*pi*24509.803921568626*t+1),timepoints)
    @test excitationField2D.excitationStrengths2 ≈ map(t->sin(2*pi*25252.52525252525*t+2),timepoints)

    excitationField3D = MPISimulations.ExcitationField(staticGradientField,sinusoidalField1,sinusoidalField2,sinusoidalField3,timepoints)
    @test excitationField3D.selectionField == staticGradientField
    @test excitationField3D.excitationField1 == sinusoidalField1.spatialFieldAmplitude
    @test excitationField3D.excitationField2 == sinusoidalField2.spatialFieldAmplitude
    @test excitationField3D.excitationField3 == sinusoidalField3.spatialFieldAmplitude
    @test excitationField3D.excitationStrengths1 ≈ map(t->sin(2*pi*24509.803921568626*t+1),timepoints)
    @test excitationField3D.excitationStrengths2 ≈ map(t->sin(2*pi*25252.52525252525*t+2),timepoints)
    @test excitationField3D.excitationStrengths3 ≈ map(t->sin(2*pi*26041.666666666668*t+3),timepoints)
end