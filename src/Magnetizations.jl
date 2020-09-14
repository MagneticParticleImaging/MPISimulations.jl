abstract type Magnetizations{T} end

# specific subtypes should implement function returning magnetization field as SVector{3,T} for given external field H, position r, and time t
@mustimplement magnetizationField(magnetization::Magnetizations{T},H::SVector{3,T},r::SVector{3,T},t::T) where T<:Real

# specific subtypes should implement function returning if magnetization is zero at position r for all times
@mustimplement isZero(magnetization::Magnetizations{T},r::SVector{3,T}) where T<:Real 

struct SeparableMagnetization{T,U<:TracerDistributions{T},V<:MagneticMoments{T}} <: Magnetizations{T}
    tracerDistribution::U
    magneticMoment::V
end

@inline isZero(separableMagnetization::SeparableMagnetization,r::SVector{3,T}) where T<:Real = isZero(separableMagnetization.tracerDistribution,r)

function magnetizationField(separableMagnetization::SeparableMagnetization{T},H::SVector{3,T},r::SVector{3,T},t::T) where T<:Real
    return tracerConcentration(separableMagnetization.tracerDistribution,r,t)*meanMagneticMoment(separableMagnetization.magneticMoment,H,t)
end