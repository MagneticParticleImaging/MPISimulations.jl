# Define abstract interface
"""
$(TYPEDEF)

Abstract supertype for tracer distributions nodes.
"""
abstract type TracerDistributions{T} end

"""
$(SIGNATURES)

Return the tracer concentration [number of particles/m³] for a given location
r [m] and time t [s].
"""
function tracerConcentration(tracerDistribution::TracerDistributions{T},r::SVector{3,T},t::T) where {T}
    error("tracerConcentration is not implemented for tracerDistribution::$(typeof(tracerDistribution)).")
end

"""
$(SIGNATURES)

Return if (`Bool`) magnetization field is zero at a given location r [m] at all
times.
"""
function isZero(tracerDistribution::TracerDistributions{T},r::SVector{3,T}) where {T}
    error("isZero is not implemented for tracerDistribution::$(typeof(tracerDistribution)).")
end

"""
$(SIGNATURES)

Return `QuadratureNodes` required to integrate over a region containing the
support of the `tracerDistribution`
"""
function QuadratureNodes(tracerDistribution::TracerDistributions{T},density) where {T}
    error("QuadratureNodes is not implemented for tracerDistribution::$(typeof(tracerDistribution)).")
end

"""
$(TYPEDEF)

Abstract supertype for static tracer distributions nodes.
"""
abstract type StaticTracerDistributions{T} <: TracerDistributions{T} end

function tracerConcentration(tracerDistribution::StaticTracerDistributions{T},position::SVector{3,T},t::T) where {T}
    if isZero(tracerDistribution,position)
        return zero(T)
    else
        return tracerDistribution.tracerConcentration
    end
end


# Define specific subtypes
## homogeneous tracer distribution within axis aligned box
struct AxisOrientedBox{T} <: StaticTracerDistributions{T}
    tracerConcentration::T
    center::SVector{3,T}
    positiveHalfLengths::SVector{3,T}
    
    AxisOrientedBox(tracerConcentration::T,center::Vector{T},positiveHalfLengths::Vector{T}) where {T} = new{T}(tracerConcentration,SVector{3,T}(center),SVector{3,T}(positiveHalfLengths))
end

@inline function isZero(tracerDistribution::AxisOrientedBox{T},r::SVector{3,T}) where {T}
    return !all(abs.(tracerDistribution.center-r).<=tracerDistribution.positiveHalfLengths)
end

function QuadratureNodes(tracerDistribution::AxisOrientedBox{T},density=10/1e-3) where {T}
    nx,ny,nz = round.(Int,2*density*tracerDistribution.positiveHalfLengths)
    nx = max(nx,1)
    ny = max(ny,1)
    nz = max(nz,1)
    n = tuple(nx,ny,nz)
    
    nodesx, weightsx = gausslegendre(nx)
    nodesx = tracerDistribution.positiveHalfLengths[1]*nodesx .+ tracerDistribution.center[1]
    weightsx *= tracerDistribution.positiveHalfLengths[1]
    nodesy, weightsy = gausslegendre(ny)
    nodesy = tracerDistribution.positiveHalfLengths[2]*nodesy .+ tracerDistribution.center[2]
    weightsy *= tracerDistribution.positiveHalfLengths[2]
    nodesz, weightsz = gausslegendre(nz)
    nodesz = tracerDistribution.positiveHalfLengths[3]*nodesz .+ tracerDistribution.center[3]
    weightsz *= tracerDistribution.positiveHalfLengths[3]   

    return GaussLegendre{T}(nodesx,weightsx,nodesy,weightsy,nodesz,weightsz,n)
end

## homogeneous tracer distribution within sphere
struct Sphere{T} <: StaticTracerDistributions{T}
    tracerConcentration::T
    center::SVector{3,T}
    radius::T
    
    Sphere(tracerConcentration::T,center::Vector{T},radius::T) where {T} = new{T}(tracerConcentration,SVector{3,T}(center),radius)
end

@inline function isZero(tracerDistribution::Sphere{T},r::SVector{3,T}) where {T}
    return norm(tracerDistribution.center-r)>tracerDistribution.radius
end

function QuadratureNodes(tracerDistribution::Sphere{T},density=10/1e-3) where {T}
    r = tracerDistribution.radius
    center = tracerDistribution.center
    lfb = center - SVector{3,T}(r,r,r)
    rrt = center + SVector{3,T}(r,r,r)
    return MidPoint(lfb,rrt,density)
end

# misc
"""
    `numberOfParticles(cmolFe=1e-3,d=25e-9)`

Most often, in experiments the iron density (mol(Fe)/L) is provided rather
than the particle concentration (number of particles per m³). This function
allows to convert the former `cmolFe` (mol(Fe)/L) to number of particles 
per m³ assuming nanoparticles with a spherical magnetite core of diameter `d`.
"""
function numberOfParticles(cmolFe=1e-3,d=25e-9)
    molarMassFe3O4 = 0.231533 # kg/mol(Fe3O4)
    densityFe3O4 = 5200 # kg/m³
    volumeNanoparticle = pi/6 * d^3 # m³
    weightNanoparticle = volumeNanoparticle .* densityFe3O4 # kg
    conversionFactor = molarMassFe3O4 / (3 * weightNanoparticle) # Number of nanoparticles per mol(Fe)/m³ (factor 3 accounts for the fact that Fe3O4 contains 3 iron atoms)
    return cmolFe * conversionFactor # Number of nanoparticles per m³
end