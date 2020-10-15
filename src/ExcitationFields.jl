abstract type ExcitationFields{T} end

# define interface
## return magnetic field strengths H as Vector{SVector{3,T}} at given location r and time t
@mustimplement excitationFieldStrength(excitation::ExcitationFields{T},r::SVector{3,T},t::T) where T<:Real


## Selection field and ND sinusoidal excitation field
struct SinusoidalExcitationField{T<:Real,N,U<:StaticMagneticFields{T}} <: ExcitationFields{T}
    selectionField::U
    excitationFields::NTuple{N,SinusoidalField{T}}
end

function SinusoidalExcitationField(selectionField::StaticMagneticFields{T}, sinusoidalFields...) where T<:Real
    excitationFields = tuple(sinusoidalFields...)
    SinusoidalExcitationField(selectionField,excitationFields)
end

function excitationFieldStrength(excitation::SinusoidalExcitationField{T},r::SVector{3,T},t::T) where T<:Real
    H = magneticFieldStrength(excitation.selectionField,r,t)
    for sinusoidalField in excitation.excitationFields
        H += magneticFieldStrength(sinusoidalField,r,t)
    end
    return H
end