module MPISimulations

using Pkg.TOML, StaticArrays, LinearAlgebra, FastGaussQuadrature

include("Misc.jl")
include("TOML.jl")

include("Quadrature.jl")
include("TracerDistributions.jl")
include("MagneticFields.jl")
include("MagneticMoments.jl")
include("Magnetizations.jl")
include("Receivers.jl")
include("ExcitationFields.jl")

end
