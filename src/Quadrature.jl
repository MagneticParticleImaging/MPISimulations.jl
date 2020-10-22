# Define abstract interface
"""
$(TYPEDEF)

Abstract supertype for quadrature nodes.
"""
abstract type QuadratureNodes{T} end

"""
$(SIGNATURES)

Implements the julia iteration interface returning (node,weight) tuples
while iterating over `quadratureNodes`.
"""
function Base.iterate(quadratureNodes::QuadratureNodes{T}) where {T}
    error("Base.iterate is not implemented for quadratureNodes::$(typeof(quadratureNodes)).")
end

function Base.iterate(quadratureNodes::QuadratureNodes{T}, state) where {T}
    error("Base.iterate is not implemented for quadratureNodes::$(typeof(quadratureNodes)).")
end


# Define specific subtypes
## gauss legendre nodes provide very good results in case, where the magnetization is continous on the integration domain
struct GaussLegendre{T} <: QuadratureNodes{T}
    nodesx::Vector{T}
    weightsx::Vector{T}
    nodesy::Vector{T}
    weightsy::Vector{T}
    nodesz::Vector{T}
    weightsz::Vector{T}
    n::Tuple{Int,Int,Int}
end

function Base.iterate(quadratureNodes::GaussLegendre{T}, state = (CartesianIndices(quadratureNodes.n),1)) where {T}
    (cartesian,counter) = state

    if counter>length(cartesian)
        return nothing
    else
        i,j,k = Tuple(cartesian[counter])
        r = SVector{3,T}(quadratureNodes.nodesx[i],quadratureNodes.nodesy[j],quadratureNodes.nodesz[k])
        w = quadratureNodes.weightsx[i]*quadratureNodes.weightsy[j]*quadratureNodes.weightsz[k]
        return ((r,w),(cartesian,counter+1))
    end
end

## midpoint nodes are the fall back option, where the magnetization is non-continous on the integration domain
struct MidPoint{T,U<:AbstractVector{T},V<:AbstractVector{T},W<:AbstractVector{T}} <: QuadratureNodes{T}
    nodesx::U
    nodesy::V
    nodesz::W
    weight::T
    n::Tuple{Int,Int,Int}
end

function Base.iterate(quadratureNodes::MidPoint{T}, state = (CartesianIndices(quadratureNodes.n),1)) where {T}
    (cartesian,counter) = state

    if counter>length(cartesian)
        return nothing
    else
        i,j,k = Tuple(cartesian[counter])
        r = SVector{3,T}(quadratureNodes.nodesx[i],quadratureNodes.nodesy[j],quadratureNodes.nodesz[k])
        return ((r,quadratureNodes.weight),(cartesian,counter+1))
    end
end

function MidPoint(lfb::SVector{3,T},rrt::SVector{3,T},density=10/1e-3) where {T}
    nx,ny,nz = round.(Int,abs.(density*(rrt-lfb)))
    nx = max(nx,1)
    ny = max(ny,1)
    nz = max(nz,1)
    n = tuple(nx,ny,nz)

    center = (rrt+lfb)/2
    halfLength = abs.((rrt-lfb)/2)    
    nodesx = range(center[1]-(1-1/nx)*halfLength[1], center[1]+(1-1/nx)*halfLength[1], length=nx)
    nodesy = range(center[2]-(1-1/ny)*halfLength[2], center[2]+(1-1/ny)*halfLength[2], length=ny)
    nodesz = range(center[3]-(1-1/nz)*halfLength[3], center[3]+(1-1/nz)*halfLength[3], length=nz)
    weight = 2*halfLength[1] * 2*halfLength[2] * 2*halfLength[3] / prod(n)

    return MidPoint{T,typeof(nodesx),typeof(nodesy),typeof(nodesz)}(nodesx,nodesy,nodesz,weight,n)
end