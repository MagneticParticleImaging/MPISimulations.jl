# Define abstract interface
"""
$(TYPEDEF)

Abstract supertype for magnetization fields.
"""
abstract type Magnetizations{T} end

"""
$(SIGNATURES)

Return the magnetization field M [T/μ₀] (μ₀ being the vacuum permeability)
as `SVector{3,T}` for a given magnetic field H [T/μ₀], location r [m] and 
time t [s].
"""
function magnetizationField(magnetization::Magnetizations{T},H::SVector{3,T},r::SVector{3,T},t::T) where {T}
    error("magnetizationField is not implemented for magnetization::$(typeof(magnetization)).")
end

"""
$(SIGNATURES)

Return if (`Bool`) magnetization field is zero at a given location r [m] at all
times.
"""
function isZero(magnetization::Magnetizations{T},r::SVector{3,T}) where {T} 
    error("isZero is not implemented for magnetization::$(typeof(magnetization)).")
end

"""
$(SIGNATURES)

Return `QuadratureNodes` required to integrate over the magnetization M.
"""
function QuadratureNodes(magnetization::Magnetizations{T},density) where {T} 
    error("QuadratureNodes is not implemented for magnetization::$(typeof(magnetization)).")
end


# Define specific subtypes
## separable magnetization
struct SeparableMagnetization{T,U<:TracerDistributions{T},V<:MagneticMoments{T}} <: Magnetizations{T}
    tracerDistribution::U
    magneticMoment::V
end

@inline isZero(magnetization::SeparableMagnetization,r::SVector{3,T}) where {T} = isZero(magnetization.tracerDistribution,r)

function magnetizationField(magnetization::SeparableMagnetization{T},H::SVector{3,T},r::SVector{3,T},t::T) where {T}
    return tracerConcentration(magnetization.tracerDistribution,r,t)*meanMagneticMoment(magnetization.magneticMoment,H,r,t)
end

function QuadratureNodes(magnetization::SeparableMagnetization{T},density=10/1e-3) where {T}
    return QuadratureNodes(magnetization.tracerDistribution,density)
end