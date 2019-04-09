module TBHamiltonians

using StaticArrays, OffsetArrays, SparseArrays
using Distributed, LinearAlgebra, Arpack
using NearestNeighbors


export Orbital, Hopping, TBHamiltonian
export sort_hoppings!, dimension, hopping_type, norbs
export add_orbital!, add_hopping!

export BlochHamiltonian, BlochHamiltonian_noF

export hamiltoniank, hamiltoniank!

export KPath, KMesh, path, mpmesh
export EigenStates, KEigenStates, KEigenAmplitudes, KBandStates, KBandAmplitudes

export bands_full, bands_sparse

export latticebasis, reciprocalbasis, spanlattice, addlattice
export cdot, toSVector

include("path_mesh_structs.jl")
include("keigen_structs.jl")


include("path_mesh_utils.jl")
include("vector_utils.jl")
include("lattice_utils.jl")


include("tbham_structs.jl")
include("tbham_utils.jl")
include("tbham_builders.jl")
include("tbham_eval.jl")


include("blochham_structs.jl")
include("blochham_utils.jl")
include("blochham_eval.jl")


include("tbham_kspectral.jl")

include("presets.jl")

end #module
