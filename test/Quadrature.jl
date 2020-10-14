using Pkg.TOML

@testset "Quadrature.jl" begin
	params = TOML.parsefile("Quadrature.toml")

    magnetization = MPISimulations.initialize(MPISimulations.SeparableMagnetization,"SeparableMagnetization1",params,Float64)
    aob = magnetization.tracerDistribution
    gl1 = MPISimulations.QuadratureNodes(magnetization)
    gl2 = MPISimulations.QuadratureNodes(aob)

    @test gl1.nodesx == gl2.nodesx
    @test gl1.weightsx == gl2.weightsx
    @test gl1.nodesy == gl2.nodesy
    @test gl1.weightsy == gl2.weightsy
    @test gl1.nodesz == gl2.nodesz
    @test gl1.weightsz == gl2.weightsz
    @test gl1.n == gl2.n
    
    for density in [1/2e-2,1/2e-3,1/2e-4,1/2e-5]
        nodesAndWeights = MPISimulations.QuadratureNodes(aob, density)
        result = 0.0
        for (r,w) in nodesAndWeights
            result += w*MPISimulations.tracerConcentration(aob,r,0.0)
        end
        @test result ≈ 4.8e-8
    end

    magnetization = MPISimulations.initialize(MPISimulations.SeparableMagnetization,"SeparableMagnetization2",params,Float64)
    sphere = magnetization.tracerDistribution
    gl1 = MPISimulations.QuadratureNodes(magnetization)
    gl2 = MPISimulations.QuadratureNodes(sphere)
    @test gl1.nodesx == gl2.nodesx
    @test gl1.nodesy == gl2.nodesy
    @test gl1.nodesz == gl2.nodesz
    @test gl1.weight == gl2.weight
    @test gl1.n == gl2.n

    for (density,rtol) in [(1/2e-2,1e0),(1/2e-3,1e0),(1/2e-4,1e-1),(1/2e-5,1e-3)]
        nodesAndWeights = MPISimulations.QuadratureNodes(sphere, density)
        result = 0.0
        for (r,w) in nodesAndWeights
            result += w*MPISimulations.tracerConcentration(sphere,r,0.0)
        end
        @test result ≈ 4/3*π*(sphere.radius)^3 rtol=rtol
    end
end