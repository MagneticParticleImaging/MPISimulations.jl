# Define abstract interface
"""
$(TYPEDEF)

Abstract supertype for magnetic moments.
"""
abstract type MagneticMoments{T} end

"""
$(SIGNATURES)

Return the mean magnetic moment [Am²] as `SVector{3,T}` for a given magnetic
field H [T/μ₀] (μ₀ being the vacuum permeability), location r [m] and time t [s].
"""
function meanMagneticMoment(magneticMoment::MagneticMoments{T},H::SVector{3,T},r::SVector{3,T},t::T) where {T}
    error("meanMagneticMoment is not implemented for magneticMoment::$(typeof(magneticMoment)).")
end


# Define specific subtypes
## Langevin model
struct Langevin{T} <: MagneticMoments{T}
    msat::T
    beta::T
end

"""
$(SIGNATURES)

Initialize `Langevin` type `MagneticMoments` with diameter `d` [m], temperature
`T` [K], relative saturation magnetization [T] `Mrel` = Msat*μ₀ (μ₀ being
the vacuum permeability) and numeric `elementType`.
"""
function Langevin(;elementType=Float64,d=25e-9,T=294,Mrel=0.6)
    k = 1.380650424e-23 #Boltzman constant
    μ₀ = 4*π*1e-7 #vacuum permeability
    M = Mrel / μ₀ #saturation magnetization of magnetic material the core is made of
    msat = M * π/6*d^3 #saturation magnetic moment of a single nanoparticle
    beta = msat /(k*T) #H measured in T/μ₀

    Langevin{elementType}(msat,beta)
end

function meanMagneticMoment(magneticMoment::Langevin{T},H::SVector{3,T},r::SVector{3,T},t::T) where {T}
    if norm(H)!=0
        x = magneticMoment.beta*norm(H)
        return magneticMoment.msat*(coth(x) - 1/x)*normalize(H)
    else
        return zero(SVector{3,T})
    end
end
