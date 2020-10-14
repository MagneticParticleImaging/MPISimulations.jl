using Pkg.TOML, StaticArrays

@testset "ExcitationFields.jl" begin
    params = TOML.parsefile("ExcitationFields.toml")
	
    Hs = MPISimulations.initialize(MPISimulations.StaticGradientField,"ExcitationField.selectionField.StaticGradientField",params,Float64)
    H0 = MPISimulations.initialize(MPISimulations.SinusoidalField,"ExcitationField.excitationField.SinusoidalField",params,Float64)
    H1 = MPISimulations.initialize(MPISimulations.SinusoidalField,"ExcitationField.excitationField1.SinusoidalField",params,Float64)
    H2 = MPISimulations.initialize(MPISimulations.SinusoidalField,"ExcitationField.excitationField2.SinusoidalField",params,Float64)
    H3 = MPISimulations.initialize(MPISimulations.SinusoidalField,"ExcitationField.excitationField3.SinusoidalField",params,Float64)
    r = randn(SVector{3,Float64})
    t = randn(Float64)
    
    excitation = MPISimulations.initialize(MPISimulations.ExcitationField1D,"ExcitationField",params,Float64)
    @test MPISimulations.excitationFieldStrength(excitation,r,t) ≈ MPISimulations.magneticFieldStrength(Hs,r,t) + MPISimulations.magneticFieldStrength(H0,r,t)

    excitation = MPISimulations.initialize(MPISimulations.ExcitationField2D,"ExcitationField",params,Float64)
    @test MPISimulations.excitationFieldStrength(excitation,r,t) ≈ MPISimulations.magneticFieldStrength(Hs,r,t) + MPISimulations.magneticFieldStrength(H1,r,t) + MPISimulations.magneticFieldStrength(H2,r,t)

    excitation = MPISimulations.initialize(MPISimulations.ExcitationField3D,"ExcitationField",params,Float64)
    @test MPISimulations.excitationFieldStrength(excitation,r,t) ≈ MPISimulations.magneticFieldStrength(Hs,r,t) + MPISimulations.magneticFieldStrength(H1,r,t) + MPISimulations.magneticFieldStrength(H2,r,t) + MPISimulations.magneticFieldStrength(H3,r,t)
end