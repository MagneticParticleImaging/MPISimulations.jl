abstract type MagneticMoments{T} end

# define interface
## return magnetic moment of a single particle as SVector{3,T} for given external field H at time t
@mustimplement meanMagneticMoment(magneticMoment::MagneticMoments{T},H::SVector{3,T},r::SVector{3,T},t::T) where T<:Real

# define specific subtypes
## Langevin Model
struct Langevin{T<:Real} <: MagneticMoments{T}
    msat::T
    beta::T
end

"""
    `initializeLangevin(S;d=25e-9,T=294,Mrel=0.6)`

Initialize Langevin type magnetic nanoparticles with diameter d (meter),
temperature T (Kelvin) and relative saturation magnetization Mrel = M*μ₀
(μ₀ being the vacuum permeability) using the keyword arguments. The underlying
numeric type can be specified via the optional argument S which can be any 
subtype of Real.
"""
function initializeLangevin(S::Type{U}=Float64;d=25e-9,T=294,Mrel=0.6) where U<:Real
    k = 1.380650424e-23 #Boltzman constant
    μ₀ = 4*π*1e-7 #vacuum permeability
    M = Mrel / μ₀ #saturation magnetization of magnetic material the core is made of
    msat = M * π/6*d^3 #saturation magnetic moment of a single nanoparticle
    beta = msat /(k*T) #H measured in T/μ₀

    Langevin{S}(msat,beta)
end

function meanMagneticMoment(magneticMoment::Langevin{T},H::SVector{3,T},r::SVector{3,T},t::T) where T<:Real
    if norm(H)!=0
        x = magneticMoment.beta*norm(H)
        return magneticMoment.msat*(coth(x) - 1/x)*normalize(H)
    else
        return zero(SVector{3,T})
    end
end
