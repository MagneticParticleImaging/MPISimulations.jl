module MPISimulations

using Pkg.TOML, StaticArrays, LinearAlgebra

include("Misc.jl")
include("TOML.jl")

include("TracerDistributions.jl")
include("MagneticFields.jl")
include("MagneticMoments.jl")
include("Magnetizations.jl")
include("Receivers.jl")

end
