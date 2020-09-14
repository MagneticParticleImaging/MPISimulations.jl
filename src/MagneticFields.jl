abstract type MagneticFields{T} end

# define interface
## return magnetic field strength H as SVector{3,T} at given location r and time t
@mustimplement magneticFieldStrength(magneticField::MagneticFields{T},r::SVector{3,T},t::T) where T<:Real


# static magnetic fields
## define abstract subtypes
abstract type StaticMagneticFields{T} <: MagneticFields{T} end

## return magnetic field strength H as SVector{3,T} at given location r
@mustimplement magneticFieldStrength(staticMagneticField::StaticMagneticFields{T},r::SVector{3,T}) where T<:Real

function magneticFieldStrength(staticMagneticField::StaticMagneticFields{T},r::SVector{3,T},t::T) where T<:Real
    return magneticFieldStrength(staticMagneticField,r)
end

## ideal gradient field
struct StaticGradientField{T<:Real} <: StaticMagneticFields{T}
    center::SVector{3,T}
    gradient::SMatrix{3,3,T}

    StaticGradientField(center::Vector{T},gradient::Matrix{T}) where T<:Real = new{T}(SVector{3,T}(center),SMatrix{3,3,T}(gradient))
end

function magneticFieldStrength(staticGradientField::StaticGradientField{T},r::SVector{3,T}) where T<:Real
    return staticGradientField.gradient*(r-staticGradientField.center)
end

## ideal homogenous field
struct StaticHomogeneousField{T<:Real} <: StaticMagneticFields{T}
    magneticFieldStrength::SVector{3,T}

    StaticHomogeneousField(magneticFieldStrength::Vector{T}) where T<:Real = new{T}(SVector{3,T}(magneticFieldStrength))
end

function magneticFieldStrength(staticHomogeneousField::StaticHomogeneousField{T},r::SVector{3,T}) where T<:Real
    return staticHomogeneousField.magneticFieldStrength
end

# dynamic magnetic fields
## ideal homogenous field with sinusoidal 1D excitation
struct SinusoidalExcitationField{T<:Real} <: MagneticFields{T}
    excitationAmplitude::StaticMagneticFields{T}
    frequency::T
    phase::T
end

function magneticFieldStrength(sinusoidalExcitationField::SinusoidalExcitationField{T},r::SVector{3,T},t::T) where T<:Real
    return magneticFieldStrength(sinusoidalExcitationField.excitationAmplitude,r)*sin(2*pi*sinusoidalExcitationField.frequency*t+sinusoidalExcitationField.phase)
end