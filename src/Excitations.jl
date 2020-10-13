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

function ExcitationField(selectionField::StaticMagneticFields{T},excitationField::SinusoidalField{T},timepoints::AbstractVector{T}) where {T<:Real}
    excitationStrengths = map(t->sin(2*pi*excitationField.frequency*t+excitationField.phase),timepoints)
    return ExcitationField1D(selectionField,excitationField.spatialFieldAmplitude,excitationStrengths)
end

function magneticExcitationFieldSequence(excitation::ExcitationField1D{T},r::SVector{3,T}) where T<:Real
    Hs = magneticFieldStrength(excitation.selectionField,r)
    H1 = magneticFieldStrength(excitation.excitationField,r)
    s1 = excitation.excitationStrengths
    fields = zeros(SVector{3,T},length(s1))

    for i in 1:length(fields)
        fields[i] = Hs + H1*s1[i]
    end
    return fields
end

## Selection field and 2D excitation field (2D MPS and 2D MPI)
struct ExcitationField2D{T<:Real,U<:StaticMagneticFields{T},V1<:StaticMagneticFields{T},V2<:StaticMagneticFields{T}} <: ExcitationFields{T}
    selectionField::U
    excitationField1::V1
    excitationField2::V2
    excitationStrengths1::Vector{T}
    excitationStrengths2::Vector{T}
end

function ExcitationField(selectionField::StaticMagneticFields{T},excitationField1::SinusoidalField{T},excitationField2::SinusoidalField{T},timepoints::AbstractVector{T}) where {T<:Real}
    excitationStrengths1 = map(t->sin(2*pi*excitationField1.frequency*t+excitationField1.phase),timepoints)
    excitationStrengths2 = map(t->sin(2*pi*excitationField2.frequency*t+excitationField2.phase),timepoints)
    return ExcitationField2D(selectionField,excitationField1.spatialFieldAmplitude,excitationField2.spatialFieldAmplitude,excitationStrengths1,excitationStrengths2)
end

function magneticExcitationFieldSequence(excitation::ExcitationField2D{T},r::SVector{3,T}) where T<:Real
    Hs = magneticFieldStrength(excitation.selectionField,r)
    H1 = magneticFieldStrength(excitation.excitationField1,r)
    H2 = magneticFieldStrength(excitation.excitationField2,r)
    s1 = excitation.excitationStrengths1
    s2 = excitation.excitationStrengths2
    fields = zeros(SVector{3,T},length(s1))

    for i in 1:length(fields)
        fields[i] = Hs + H1*s1[i] + H2*s2[i]
    end
    return fields
end

## Selection field and 3D excitation field (3D MPS and 3D MPI)
struct ExcitationField3D{T<:Real,U<:StaticMagneticFields{T},V1<:StaticMagneticFields{T},V2<:StaticMagneticFields{T},V3<:StaticMagneticFields{T}} <: ExcitationFields{T}
    selectionField::U
    excitationField1::V1
    excitationField2::V2
    excitationField3::V3
    excitationStrengths1::Vector{T}
    excitationStrengths2::Vector{T}
    excitationStrengths3::Vector{T}
end

function ExcitationField(selectionField::StaticMagneticFields{T},excitationField1::SinusoidalField{T},excitationField2::SinusoidalField{T},excitationField3::SinusoidalField{T},timepoints::AbstractVector{T}) where {T<:Real}
    excitationStrengths1 = map(t->sin(2*pi*excitationField1.frequency*t+excitationField1.phase),timepoints)
    excitationStrengths2 = map(t->sin(2*pi*excitationField2.frequency*t+excitationField2.phase),timepoints)
    excitationStrengths3 = map(t->sin(2*pi*excitationField3.frequency*t+excitationField3.phase),timepoints)
    return ExcitationField3D(selectionField,excitationField1.spatialFieldAmplitude,excitationField2.spatialFieldAmplitude,excitationField3.spatialFieldAmplitude,excitationStrengths1,excitationStrengths2,excitationStrengths3)
end

function magneticExcitationFieldSequence(excitation::ExcitationField3D{T},r::SVector{3,T}) where T<:Real
    Hs = magneticFieldStrength(excitation.selectionField,r)
    H1 = magneticFieldStrength(excitation.excitationField1,r)
    H2 = magneticFieldStrength(excitation.excitationField2,r)
    H3 = magneticFieldStrength(excitation.excitationField3,r)
    s1 = excitation.excitationStrengths1
    s2 = excitation.excitationStrengths2
    s3 = excitation.excitationStrengths3
    fields = zeros(SVector{3,T},length(s1))

    for i in 1:length(fields)
        fields[i] = Hs + H1*s1[i] + H2*s2[i] + H3*s3[i]
    end
    return fields
end