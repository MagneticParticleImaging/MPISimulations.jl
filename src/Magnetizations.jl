abstract type Magnetizations{T} end

# specific subtypes should implement function returning magnetization field as SVector{3,T} for given external field H, position r, and time t
@mustimplement magnetizationField(magnetization::Magnetizations{T},H::SVector{3,T},r::SVector{3,T},t::T) where T<:Real

# specific subtypes should implement function returning if magnetization is zero at position r for all times
@mustimplement isZero(magnetization::Magnetizations{T},r::SVector{3,T}) where T<:Real 

# specific subtypes should implement function returning quadrature nodes
@mustimplement quadratureNodes(magnetization::Magnetizations{T},density) where T<:Real

# separable magnetization
struct SeparableMagnetization{T,U<:TracerDistributions{T},V<:MagneticMoments{T}} <: Magnetizations{T}
    tracerDistribution::U
    magneticMoment::V
end

@inline isZero(magnetization::SeparableMagnetization,r::SVector{3,T}) where T<:Real = isZero(magnetization.tracerDistribution,r)

function magnetizationField(magnetization::SeparableMagnetization{T},H::SVector{3,T},r::SVector{3,T},t::T) where T<:Real
    return tracerConcentration(magnetization.tracerDistribution,r,t)*meanMagneticMoment(magnetization.magneticMoment,H,r,t)
end

function quadratureNodes(magnetization::SeparableMagnetization{T},density=10/1e-3) where T<:Real
    return quadratureNodes(magnetization.tracerDistribution,density)
end