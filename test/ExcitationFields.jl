using Pkg.TOML, StaticArrays

@testset "ExcitationFields.jl" begin
    params = TOML.parsefile("ExcitationFields.toml")
	
    Hs = MPISimulations.initialize(MPISimulations.StaticGradientField,"StaticGradientField",params,Float64)
    H1 = MPISimulations.initialize(MPISimulations.SinusoidalField,"SinusoidalField1",params,Float64)
    H2 = MPISimulations.initialize(MPISimulations.SinusoidalField,"SinusoidalField2",params,Float64)
    H3 = MPISimulations.initialize(MPISimulations.SinusoidalField,"SinusoidalField3",params,Float64)
    r = randn(SVector{3,Float64})
    t = randn(Float64)
    
    excitation = MPISimulations.ExcitationField(Hs,H1)
    @test MPISimulations.magneticFieldStrength(excitation,r,t) ≈ MPISimulations.magneticFieldStrength(Hs,r,t) + MPISimulations.magneticFieldStrength(H1,r,t)

    excitation = MPISimulations.ExcitationField(Hs,H1,H2)
    @test MPISimulations.magneticFieldStrength(excitation,r,t) ≈ MPISimulations.magneticFieldStrength(Hs,r,t) + MPISimulations.magneticFieldStrength(H1,r,t) + MPISimulations.magneticFieldStrength(H2,r,t)

    excitation = MPISimulations.ExcitationField(Hs,H1,H2,H3)
    @test MPISimulations.magneticFieldStrength(excitation,r,t) ≈ MPISimulations.magneticFieldStrength(Hs,r,t) + MPISimulations.magneticFieldStrength(H1,r,t) + MPISimulations.magneticFieldStrength(H2,r,t) + MPISimulations.magneticFieldStrength(H3,r,t)
end