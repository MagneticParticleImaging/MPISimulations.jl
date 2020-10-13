abstract type ExcitationFields{T} end

# define interface
## return magnetic field strengths H as Vector{SVector{3,T}} at given location r
@mustimplement magneticExcitationFieldSequence(excitation::ExcitationFields{T},r::SVector{3,T}) where T<:Real


## Selection field and 1D excitation field (1D MPS and 1D MPI)
struct ExcitationField1D{T<:Real,U<:StaticMagneticFields{T},V<:StaticMagneticFields{T}} <: ExcitationFields{T}
    selectionField::U
    excitationField::V
    excitationStrengths::Vector{T}
end

function ExcitationField1D(selectionField::StaticMagneticFields{T},excitationField::SinusoidalField{T},timepoints::AbstractVector{T}) where {T<:Real}
    excitationStrengths = map(t->sin(2*pi*excitationField.frequency*t+excitationField.phase),timepoints)
    return ExcitationField1D(selectionField,excitationField.spatialFieldAmplitude,excitationStrengths)
end

function magneticExcitationFieldSequence(excitation::ExcitationField1D{T},r::SVector{3,T}) where T<:Real
    H₀ = magneticFieldStrength(excitation.selectionField,r)
    Hₓ = magneticFieldStrength(excitation.excitationField,r)
    s = excitation.excitationStrengths
    fields = zeros(SVector{3,T},length(s))

    for i in 1:length(fields)
        fields[i] = H₀ + Hₓ*s[i]
    end
    return fields
end