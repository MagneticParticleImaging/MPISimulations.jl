# Define abstract interface
"""
$(TYPEDEF)

Abstract supertype for excitation fields.
"""
abstract type ExcitationFields{T,N} end

"""
$(SIGNATURES)

Initialize an excitation field given a `selectionField` and multiple 
`excitationFields`.
"""
function ExcitationField(selectionField::MagneticFields{T}, excitationFields::U) where {T,U<:Tuple{Vararg{W} where W<:MagneticFields{T}}}
    error("ExcitationField is not implemented for selectionField::$(typeof(selectionField)) and excitationFields::$(typeof(excitationFields)).")
end

function ExcitationField(selectionField::MagneticFields, excitationFields...)
    return ExcitationField(selectionField,tuple(excitationFields...))
end

"""
$(SIGNATURES)

Return the magnetic field strength H [T/μ₀] (μ₀ being the vacuum permeability)
as `SVector{3,T}` for a given location r [m] and time t [s].
"""
function magneticFieldStrength(excitation::ExcitationFields{T,N},r::SVector{3,T},t::T) where {T,N}
    error("magneticFieldStrength is not implemented for $(typeof(excitation)).")
end


# Define specific subtypes
## Selection field and N-dimensional sinusoidal excitation field
struct SinusoidalExcitationField{T,N,U<:StaticMagneticFields{T}} <: ExcitationFields{T,N}
    selectionField::U
    sinusoidalFields::NTuple{N,SinusoidalField{T}}
end

function ExcitationField(selectionField::U, magneticFields::NTuple{N,SinusoidalField{T}}) where {T,N,U<:StaticMagneticFields{T}}
    return SinusoidalExcitationField{T,N,U}(selectionField,magneticFields)
end

function magneticFieldStrength(excitation::SinusoidalExcitationField{T,N},r::SVector{3,T},t::T) where {T,N}
    H = magneticFieldStrength(excitation.selectionField,r,t)
    @inbounds for i =1:N
        H += magneticFieldStrength(excitation.sinusoidalFields[i],r,t)
    end
    return H
end