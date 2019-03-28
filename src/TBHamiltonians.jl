module TBHamiltonians

using StaticArrays, OffsetArrays, LinearAlgebra, SparseArrays

export Orbital, Hopping, TBHamiltonian
export sort_hoppings!, dimension, hopping_type, norbs
export add_orbital!, add_hopping!

export BlochHamiltonian

export hamiltoniank, hamiltoniank!

export KPath, KMesh

export cdot, toSVector



include("tbham_structs.jl")
include("tbham_utils.jl")
include("tbham_builders.jl")
include("vector_utils.jl")
include("evalhamiltonian.jl")

include("path_mesh.jl")
include("path_mesh_utils.jl")

include("blochham_structs.jl")
include("blochham_builders.jl")

end #module
