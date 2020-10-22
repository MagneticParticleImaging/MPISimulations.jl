using Pkg.TOML, StaticArrays

@testset "Receivers.jl" begin
    r = rand(SVector{3,Float64})
    t = rand(Float64)

    params = TOML.parsefile("Receivers.toml")
    
    xSensitivity = MPISimulations.initialize(MPISimulations.StaticHomogeneousField,"Receiver.xSensitivity",params,Float64)
    ySensitivity = MPISimulations.initialize(MPISimulations.StaticHomogeneousField,"Receiver.ySensitivity",params,Float64)
    zSensitivity = MPISimulations.initialize(MPISimulations.StaticHomogeneousField,"Receiver.zSensitivity",params,Float64)


    samplingRate = params["Receiver"]["samplingRate"]
    
    receiver1D = MPISimulations.ReceiveCoils(samplingRate,xSensitivity)
    @test MPISimulations.receiveCoilSensitivities(receiver1D,r,t) == SMatrix{1,3}(0.1,0.0,0.0)

    receiver2D = MPISimulations.ReceiveCoils(samplingRate,xSensitivity,ySensitivity)
    @test MPISimulations.receiveCoilSensitivities(receiver2D,r,t) == SMatrix{2,3}(0.1,0.0,0.0,0.11,0.0,0.0)

    receiver3D = MPISimulations.ReceiveCoils(samplingRate,xSensitivity,ySensitivity,zSensitivity)
    @test MPISimulations.receiveCoilSensitivities(receiver3D,r,t) == SMatrix{3,3}(0.1,0.0,0.0,0.0,0.11,0.0,0.0,0.0,0.12)
end