using Pkg.TOML, StaticArrays

@testset "Receivers.jl" begin
    r = rand(SVector{3,Float64})

    params = TOML.parsefile("Receivers.toml")
	
    receiver1D = MPISimulations.initialize(MPISimulations.Receiver1D,"Receiver",params,Float64)
    @test MPISimulations.receiverCoilSensitivities(receiver1D,r) == SMatrix{1,3}(0.12,0.12,0.12)

    receiver2D = MPISimulations.initialize(MPISimulations.Receiver2D,"Receiver",params,Float64)
    @test MPISimulations.receiverCoilSensitivities(receiver2D,r) == SMatrix{2,3}(0.12,0.13,0.12,0.13,0.12,0.13)

    receiver3D = MPISimulations.initialize(MPISimulations.Receiver3D,"Receiver",params,Float64)
    @test MPISimulations.receiverCoilSensitivities(receiver3D,r) == SMatrix{3,3}(0.12,0.13,0.14,0.12,0.13,0.14,0.12,0.13,0.14)
end