abstract type ExcitationFields{T} end

# define interface
## return magnetic field strengths H as Vector{SVector{3,T}} at given location r and time t
@mustimplement excitationFieldStrength(excitation::ExcitationFields{T},r::SVector{3,T},t::T) where T<:Real


## Selection field and 1D excitation field (1D MPS and 1D MPI)
struct ExcitationField1D{T<:Real,U<:StaticMagneticFields{T},V<:MagneticFields{T}} <: ExcitationFields{T}
    selectionField::U
    excitationField::V
end

function excitationFieldStrength(excitation::ExcitationField1D{T},r::SVector{3,T},t::T) where T<:Real
    Hs = magneticFieldStrength(excitation.selectionField,r,t)
    H1 = magneticFieldStrength(excitation.excitationField,r,t)
    return Hs + H1
end

## Selection field and 2D excitation field (2D MPS and 2D MPI)
struct ExcitationField2D{T<:Real,U<:StaticMagneticFields{T},V<:MagneticFields{T},W<:MagneticFields{T}} <: ExcitationFields{T}
    selectionField::U
    excitationField1::V
    excitationField2::W
end

function excitationFieldStrength(excitation::ExcitationField2D{T},r::SVector{3,T},t::T) where T<:Real
    Hs = magneticFieldStrength(excitation.selectionField,r,t)
    H1 = magneticFieldStrength(excitation.excitationField1,r,t)
    H2 = magneticFieldStrength(excitation.excitationField2,r,t)
    return Hs + H1 + H2
end

## Selection field and 3D excitation field (3D MPS and 3D MPI)
struct ExcitationField3D{T<:Real,U<:StaticMagneticFields{T},V<:MagneticFields{T},W<:MagneticFields{T},X<:MagneticFields{T}} <: ExcitationFields{T}
    selectionField::U
    excitationField1::V
    excitationField2::W
    excitationField3::X
end

function excitationFieldStrength(excitation::ExcitationField3D{T},r::SVector{3,T},t::T) where T<:Real
    Hs = magneticFieldStrength(excitation.selectionField,r,t)
    H1 = magneticFieldStrength(excitation.excitationField1,r,t)
    H2 = magneticFieldStrength(excitation.excitationField2,r,t)
    H3 = magneticFieldStrength(excitation.excitationField3,r,t)
    return Hs + H1 + H2 + H3
end