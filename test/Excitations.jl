using Pkg.TOML, StaticArrays

@testset "Excitations.jl" begin
    params = TOML.parsefile("Excitations.toml")
	
    excitationField1D = MPISimulations.initialize(MPISimulations.ExcitationField1D,"ExcitationField1D",params,Float64)
    r = SVector{3,Float64}(0.0,0.0,0.1)
    @test MPISimulations.magneticExcitationFieldSequence(excitationField1D,r) == SVector{3,Float64}[SVector{3,Float64}(0.12, 0.0, -0.2),SVector{3,Float64}(0.0, 0.0, -0.2),SVector{3,Float64}(-0.12, 0.0, -0.2)]
    r = zero(SVector{3,Float64})
    @test MPISimulations.magneticExcitationFieldSequence(excitationField1D,r) == SVector{3,Float64}[SVector{3,Float64}(0.12, 0.0, 0.0),SVector{3,Float64}(0.0, 0.0, 0.0),SVector{3,Float64}(-0.12, 0.0, 0.0)]

    # create excitation from static and sinusoidal field
    staticGradientField = MPISimulations.initialize(MPISimulations.StaticGradientField,"StaticGradientField",params,Float64)
    sinusoidalField1 = MPISimulations.initialize(MPISimulations.SinusoidalField,"SinusoidalField1",params,Float64)
    sinusoidalField2 = MPISimulations.initialize(MPISimulations.SinusoidalField,"SinusoidalField2",params,Float64)
    sinusoidalField3 = MPISimulations.initialize(MPISimulations.SinusoidalField,"SinusoidalField3",params,Float64)
    timepoints = range(0, 1, length=5)
    Hs = MPISimulations.magneticFieldStrength(staticGradientField,r)
    H1 = MPISimulations.magneticFieldStrength(sinusoidalField1.spatialFieldAmplitude,r)
    H2 = MPISimulations.magneticFieldStrength(sinusoidalField2.spatialFieldAmplitude,r)
    H3 = MPISimulations.magneticFieldStrength(sinusoidalField3.spatialFieldAmplitude,r)
    f1 = 24509.803921568626
    f2 = 25252.52525252525
    f3 = 26041.666666666668
    phase1 = 1.0
    phase2 = 2.0
    phase3 = 3.0

    excitationField1D = MPISimulations.ExcitationField(staticGradientField,sinusoidalField1,timepoints)
    @test excitationField1D.selectionField == staticGradientField
    @test excitationField1D.excitationField == sinusoidalField1.spatialFieldAmplitude
    @test excitationField1D.excitationStrengths ≈ map(t->sin(2*pi*f1*t+1phase1),timepoints)
    @test MPISimulations.magneticExcitationFieldSequence(excitationField1D,r) == map(t->Hs+sin(2*pi*f1*t+phase1)*H1,timepoints)

    excitationField2D = MPISimulations.ExcitationField(staticGradientField,sinusoidalField1,sinusoidalField2,timepoints)
    @test excitationField2D.selectionField == staticGradientField
    @test excitationField2D.excitationField1 == sinusoidalField1.spatialFieldAmplitude
    @test excitationField2D.excitationField2 == sinusoidalField2.spatialFieldAmplitude
    @test excitationField2D.excitationStrengths1 ≈ map(t->sin(2*pi*f1*t+phase1),timepoints)
    @test excitationField2D.excitationStrengths2 ≈ map(t->sin(2*pi*f2*t+phase2),timepoints)
    @test MPISimulations.magneticExcitationFieldSequence(excitationField2D,r) == map(t->Hs+sin(2*pi*f1*t+phase1)*H1+sin(2*pi*f2*t+phase2)*H2,timepoints)

    excitationField3D = MPISimulations.ExcitationField(staticGradientField,sinusoidalField1,sinusoidalField2,sinusoidalField3,timepoints)
    @test excitationField3D.selectionField == staticGradientField
    @test excitationField3D.excitationField1 == sinusoidalField1.spatialFieldAmplitude
    @test excitationField3D.excitationField2 == sinusoidalField2.spatialFieldAmplitude
    @test excitationField3D.excitationField3 == sinusoidalField3.spatialFieldAmplitude
    @test excitationField3D.excitationStrengths1 ≈ map(t->sin(2*pi*f1*t+phase1),timepoints)
    @test excitationField3D.excitationStrengths2 ≈ map(t->sin(2*pi*f2*t+phase2),timepoints)
    @test excitationField3D.excitationStrengths3 ≈ map(t->sin(2*pi*f3*t+phase3),timepoints)
    @test MPISimulations.magneticExcitationFieldSequence(excitationField3D,r) == map(t->Hs+sin(2*pi*f1*t+phase1)*H1+sin(2*pi*f2*t+phase2)*H2+sin(2*pi*f3*t+phase3)*H3,timepoints)
end