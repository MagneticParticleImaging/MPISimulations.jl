abstract type TracerDistributions{T} end

# define interface
## return tracer concentration as as type T fro given location r and time t
@mustimplement tracerConcentration(tracerDistribution::TracerDistributions{T},r::SVector{3,T},t::T) where T<:Real

# define specific subtypes
## static tracer distributions
abstract type StaticTracerDistributions{T} <: TracerDistributions{T} end

## define interface for static tracer distributions
@mustimplement isInside(staticTracerDistributions::StaticTracerDistributions{T},r::SVector{3,T}) where T<:Real

function tracerConcentration(tracerDistribution::StaticTracerDistributions{T},position::SVector{3,T},t::T) where T<:Real
    if isInside(tracerDistribution,position)
        return tracerDistribution.tracerConcentration
    else
        return zero(T)
    end
end

### homogeneous tracer distribution within axis aligned box
struct AxisOrientedBox{T<:Real} <: StaticTracerDistributions{T}
    tracerConcentration::T
    center::SVector{3,T}
    positiveHalfLengths::SVector{3,T}
    
    AxisOrientedBox(tracerConcentration::T,center::Vector{T},positiveHalfLengths::Vector{T}) where T<:Real = new{T}(tracerConcentration,SVector{3,T}(center),SVector{3,T}(positiveHalfLengths))
end

@inline function isInside(axisOrientedBox::AxisOrientedBox{T},r::SVector{3,T}) where T<:Real
    return all(abs.(axisOrientedBox.center-r).<=axisOrientedBox.positiveHalfLengths)
end

### homogeneous tracer distribution within sphere
struct Sphere{T<:Real} <: StaticTracerDistributions{T}
    tracerConcentration::T
    center::SVector{3,T}
    radius::T
    
    Sphere(tracerConcentration::T,center::Vector{T},radius::T) where T<:Real = new{T}(tracerConcentration,SVector{3,T}(center),radius)
end

@inline function isInside(sphere::Sphere{T},r::SVector{3,T}) where T<:Real
    return norm(sphere.center-r)<=sphere.radius
end
