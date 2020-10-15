abstract type MagneticFields{T} end

# define interface
## return magnetic field strength H as SVector{3,T} at given location r and time t
@mustimplement magneticFieldStrength(magneticField::MagneticFields{T},r::SVector{3,T},t::T) where T<:Real


# static magnetic fields
## define abstract subtypes
abstract type StaticMagneticFields{T} <: MagneticFields{T} end

## return magnetic field strength H as SVector{3,T} at given location r
@mustimplement magneticFieldStrength(magneticField::StaticMagneticFields{T},r::SVector{3,T}) where T<:Real

function magneticFieldStrength(magneticField::StaticMagneticFields{T},r::SVector{3,T},t::T) where T<:Real
    return magneticFieldStrength(magneticField,r)
end

## ideal gradient field
struct StaticGradientField{T<:Real} <: StaticMagneticFields{T}
    center::SVector{3,T}
    gradient::SMatrix{3,3,T}

    StaticGradientField(center::AbstractVector{T},gradient::AbstractMatrix{T}) where T<:Real = new{T}(SVector{3,T}(center),SMatrix{3,3,T}(gradient))
end

function magneticFieldStrength(magneticField::StaticGradientField{T},r::SVector{3,T}) where T<:Real
    return magneticField.gradient*(r-magneticField.center)
end

## ideal homogenous field
struct StaticHomogeneousField{T<:Real} <: StaticMagneticFields{T}
    magneticFieldStrength::SVector{3,T}

    StaticHomogeneousField(magneticFieldStrength::Vector{T}) where T<:Real = new{T}(SVector{3,T}(magneticFieldStrength))
end

function magneticFieldStrength(magneticField::StaticHomogeneousField{T},r::SVector{3,T}) where T<:Real
    return magneticField.magneticFieldStrength
end

# dynamic magnetic fields
## ideal homogenous field with sinusoidal 1D excitation
struct SinusoidalField{T<:Real} <: MagneticFields{T}
    spatialFieldProfile::StaticMagneticFields{T}
    baseFrequency::T
    divider::Int
    doubleFrequency::T
    relativePhase::T # = phase/π
end

function SinusoidalField(spatialFieldProfile::StaticMagneticFields{T},baseFrequency::Real,divider::Real,relativePhase::Real=0.0) where T<:Real
    doubleFrequency = 2//divider * baseFrequency
    return SinusoidalField{T}(spatialFieldProfile,T(baseFrequency),T(divider),T(doubleFrequency),T(relativePhase))
end

function magneticFieldStrength(magneticField::SinusoidalField{T},r::SVector{3,T},t::T) where T<:Real
    return magneticFieldStrength(magneticField.spatialFieldProfile,r)*sinpi(magneticField.doubleFrequency*t+magneticField.relativePhase)
end