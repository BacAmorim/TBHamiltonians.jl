module TBHamiltonians

using StaticArrays, OffsetArrays, LinearAlgebra

export Orbital, Hopping, TBHamiltonian
export dimension,  hopping_type, norbs
export add_orbital!, add_hopping!

include("tbham_structs.jl")
include("tbham_utils.jl")
include("tbham_builders.jl")

end #module
