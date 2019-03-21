module TBHamiltonians

using StaticArrays, OffsetArrays, LinearAlgebra, SparseArrays

export Orbital, Hopping, TBHamiltonian
export sort_hoppings!, dimension, hopping_type, norbs
export cdot
export add_orbital!, add_hopping!
export hamiltoniank, hamiltoniank!

include("tbham_structs.jl")
include("tbham_utils.jl")
include("tbham_builders.jl")
include("vectorutils.jl")
include("evalhamiltonian.jl")

end #module
