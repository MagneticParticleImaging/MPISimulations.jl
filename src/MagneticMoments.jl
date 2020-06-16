abstract type MagneticMoments{T} end

# define interface
## return magnetic moment as SVector{3,T} for given external field H at time t
@mustimplement meanMagneticMoment(magneticMoment::MagneticMoments{T},H::SVector{3,T},t::T) where T<:Real

# define specific subtypes
## Langevin Model
struct Langevin{T<:Real} <: MagneticMoments{T}
    msat::T
    beta::T
end

function meanMagneticMoment(langevin::Langevin{T},H::SVector{3,T},t::T) where T<:Real
    if norm(H)!=0
        x = langevin.beta*norm(H)
        return langevin.msat*(coth(x) - 1/x)*normalize(H)
    else
        return zero(SVector{3,T})
    end
end
