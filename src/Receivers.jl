abstract type Receivers{T} end

# define interface
## return receive coil sensitivity (3 components each) of all L receive
## channels at given location r as Lx3 matrix.
@mustimplement receiverCoilSensitivities(receiver::Receivers{T},r::SVector{3,T}) where T<:Real

# define specific subtypes
## 1D receiver
struct Receiver1D{T,U<:StaticMagneticFields{T}} <: Receivers{T}
    samplingRate::Real
    xSensitivity::U
end

function receiverCoilSensitivities(receiver::Receiver1D{T},r::SVector{3,T}) where T<:Real
    px = magneticFieldStrength(receiver.xSensitivity,r)

    return SMatrix{1,3}(px)
end

## 2D receiver
struct Receiver2D{T,U<:StaticMagneticFields{T},V<:StaticMagneticFields{T}} <: Receivers{T}
    samplingRate::Real
    xSensitivity::U
    ySensitivity::V
end

function receiverCoilSensitivities(receiver::Receiver2D{T},r::SVector{3,T}) where T<:Real
    px = magneticFieldStrength(receiver.xSensitivity,r)
    py = magneticFieldStrength(receiver.ySensitivity,r)

    return transpose(SMatrix{3,2}(px...,py...))
end

## 3D receiver
struct Receiver3D{T,U<:StaticMagneticFields{T},V<:StaticMagneticFields{T},W<:StaticMagneticFields{T}} <: Receivers{T}
    samplingRate::Real
    xSensitivity::U
    ySensitivity::V
    zSensitivity::W
end

function receiverCoilSensitivities(receiver::Receiver3D{T},r::SVector{3,T}) where T<:Real
    px = magneticFieldStrength(receiver.xSensitivity,r)
    py = magneticFieldStrength(receiver.ySensitivity,r)
    pz = magneticFieldStrength(receiver.zSensitivity,r)

    return transpose(SMatrix{3,3}(px...,py...,pz...))
end