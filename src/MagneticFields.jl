# Define abstract interface
"""
$(TYPEDEF)

Abstract supertype for magnetic fields.
"""
abstract type MagneticFields{T} end

"""
$(SIGNATURES)

Return the magnetic field strength H [T/μ₀] (μ₀ being the vacuum permeability)
as `SVector{3,T}` for a given location r [m] and time t [s].
"""
function magneticFieldStrength(magneticField::MagneticFields{T},r::SVector{3,T},t::T) where {T}
    error("magneticFieldStrength is not implemented for magneticField::$(typeof(magneticField)).")
end

"""
$(TYPEDEF)

Abstract supertype for static magnetic fields.
"""
abstract type StaticMagneticFields{T} <: MagneticFields{T} end

"""
$(SIGNATURES)

Return the magnetic field strength H [T/μ₀] (μ₀ being the vacuum permeability)
as `SVector{3,T}` for a given location r [m].
"""
function magneticFieldStrength(magneticField::StaticMagneticFields{T},r::SVector{3,T}) where {T}
    error("magneticFieldStrength is not implemented for magneticField::$(typeof(magneticField)).")
end

function magneticFieldStrength(magneticField::StaticMagneticFields{T},r::SVector{3,T},t::T) where {T}
    return magneticFieldStrength(magneticField,r)
end


# Define specific subtypes
## ideal gradient field
struct StaticGradientField{T} <: StaticMagneticFields{T}
    center::SVector{3,T}
    gradient::SMatrix{3,3,T}

    StaticGradientField(center::AbstractVector{T},gradient::AbstractMatrix{T}) where {T} = new{T}(SVector{3,T}(center),SMatrix{3,3,T}(gradient))
end

function magneticFieldStrength(magneticField::StaticGradientField{T},r::SVector{3,T}) where {T}
    return magneticField.gradient*(r-magneticField.center)
end

## ideal homogenous field
struct StaticHomogeneousField{T} <: StaticMagneticFields{T}
    magneticFieldStrength::SVector{3,T}

    StaticHomogeneousField(magneticFieldStrength::Vector{T}) where {T} = new{T}(SVector{3,T}(magneticFieldStrength))
end

function magneticFieldStrength(magneticField::StaticHomogeneousField{T},r::SVector{3,T}) where {T}
    return magneticField.magneticFieldStrength
end

## ideal homogenous field with sinusoidal 1D excitation
struct SinusoidalField{T} <: MagneticFields{T}
    spatialFieldProfile::StaticMagneticFields{T}
    baseFrequency::T
    divider::Int
    doubleFrequency::T
    relativePhase::T # = phase/π
end

function SinusoidalField(spatialFieldProfile::StaticMagneticFields{T},baseFrequency::Real,divider::Real,relativePhase::Real=0.0) where {T}
    doubleFrequency = 2//divider * baseFrequency
    return SinusoidalField{T}(spatialFieldProfile,T(baseFrequency),T(divider),T(doubleFrequency),T(relativePhase))
end

function magneticFieldStrength(magneticField::SinusoidalField{T},r::SVector{3,T},t::T) where {T}
    return magneticFieldStrength(magneticField.spatialFieldProfile,r)*sinpi(magneticField.doubleFrequency*t+magneticField.relativePhase)
end