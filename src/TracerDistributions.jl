abstract type TracerDistributions{T} end

# define interface
## specific subtypes should implement function returning number of nanoparticles per m³ as as type T for given location r and time t
@mustimplement tracerConcentration(tracerDistribution::TracerDistributions{T},r::SVector{3,T},t::T) where T<:Real

## specific subtypes should implement function returning if tracer distribution is zero at position r for all times
@mustimplement isZero(tracerDistributions::TracerDistributions{T},r::SVector{3,T}) where T<:Real

# define abstract subtypes
## static tracer distributions
abstract type StaticTracerDistributions{T} <: TracerDistributions{T} end

function tracerConcentration(tracerDistribution::StaticTracerDistributions{T},position::SVector{3,T},t::T) where T<:Real
    if isZero(tracerDistribution,position)
        return zero(T)
    else
        return tracerDistribution.tracerConcentration
    end
end

### homogeneous tracer distribution within axis aligned box
struct AxisOrientedBox{T<:Real} <: StaticTracerDistributions{T}
    tracerConcentration::T
    center::SVector{3,T}
    positiveHalfLengths::SVector{3,T}
    
    AxisOrientedBox(tracerConcentration::T,center::Vector{T},positiveHalfLengths::Vector{T}) where T<:Real = new{T}(tracerConcentration,SVector{3,T}(center),SVector{3,T}(positiveHalfLengths))
end

@inline function isZero(tracerDistribution::AxisOrientedBox{T},r::SVector{3,T}) where T<:Real
    return !all(abs.(tracerDistribution.center-r).<=tracerDistribution.positiveHalfLengths)
end

### homogeneous tracer distribution within sphere
struct Sphere{T<:Real} <: StaticTracerDistributions{T}
    tracerConcentration::T
    center::SVector{3,T}
    radius::T
    
    Sphere(tracerConcentration::T,center::Vector{T},radius::T) where T<:Real = new{T}(tracerConcentration,SVector{3,T}(center),radius)
end

@inline function isZero(tracerDistribution::Sphere{T},r::SVector{3,T}) where T<:Real
    return norm(tracerDistribution.center-r)>tracerDistribution.radius
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