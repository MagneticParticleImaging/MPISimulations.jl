# Define abstract interface
"""
$(TYPEDEF)

Abstract supertype for receive coils.
"""
abstract type ReceiveCoils{T,N} end

"""
$(SIGNATURES)

Initialize an receive coil setup with `samplingRate` [Hz] and multiple receive
coils specified by their coil sensitivity [ToDo].
`excitationFields`.
"""
function ReceiveCoils(samplingRate::Real, sensitivities::U) where {T,U<:Tuple{Vararg{W} where W<:MagneticFields{T}}}
    error("ExcitationField is not implemented for selectionField::$(typeof(selectionField)) and excitationFields::$(typeof(excitationFields)).")
end

function ReceiveCoils(samplingRate::Real, sensitivities...)
    return ReceiveCoils(samplingRate,tuple(sensitivities...))
end

"""
$(SIGNATURES)

Return a receive coil sensitivities [??] as SMatrix{N,3,T} for a given location
r [m] and time t [s].
"""
function receiveCoilSensitivities(receiveCoil::ReceiveCoils{T,N},r::SVector{3,T},t::T) where {T,N}
    error("receiveCoilSensitivities is not implemented for receiveCoil::$(typeof(receiveCoil)).")
end


# define specific subtypes
## N-dimensional receive coil setup all running with the same syncronized clock
struct StaticReceiveCoils{T,N,V<:Tuple{Vararg{W} where W<:StaticMagneticFields{T}},U} <: ReceiveCoils{T,N}
    samplingRate::U
    sensitivities::V
end

function receiveCoilSensitivities(receiveCoil::StaticReceiveCoils{T,N},r::SVector{3,T},t::T) where {T,N}
    m = MMatrix{N,3,T}(undef)
    @inbounds for i =1:N
        m[i,:] = magneticFieldStrength(receiveCoil.sensitivities[i],r,t)
    end
    return SMatrix{N,3,T}(m)
end

function ReceiveCoils(samplingRate::Real,sensitivities::V) where {T,V<:Tuple{Vararg{W} where W<:StaticMagneticFields{T}}}
    N = length(sensitivities)
    return StaticReceiveCoils{T,N,V,typeof(samplingRate)}(samplingRate,sensitivities)
end