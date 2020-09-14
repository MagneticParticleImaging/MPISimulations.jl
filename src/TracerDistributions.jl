abstract type TracerDistributions{T} end

# define interface
## specific subtypes should implement function returning tracer concentration as as type T for given location r and time t
@mustimplement tracerConcentration(tracerDistribution::TracerDistributions{T},r::SVector{3,T},t::T) where T<:Real

## specific subtypes should implement function returning if tracer distribution is zero at position r for all times
@mustimplement isZero(tracerDistributions::TracerDistributions{T},r::SVector{3,T}) where T<:Real

# define specific subtypes
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

@inline function isZero(axisOrientedBox::AxisOrientedBox{T},r::SVector{3,T}) where T<:Real
    return !all(abs.(axisOrientedBox.center-r).<=axisOrientedBox.positiveHalfLengths)
end

### homogeneous tracer distribution within sphere
struct Sphere{T<:Real} <: StaticTracerDistributions{T}
    tracerConcentration::T
    center::SVector{3,T}
    radius::T
    
    Sphere(tracerConcentration::T,center::Vector{T},radius::T) where T<:Real = new{T}(tracerConcentration,SVector{3,T}(center),radius)
end

@inline function isZero(sphere::Sphere{T},r::SVector{3,T}) where T<:Real
    return norm(sphere.center-r)>sphere.radius
end
