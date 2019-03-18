module TBHamiltonians

using StaticArrays, LinearAlgebra

export Orbital, Hopping, TBHamiltonian
export dimension,  hopping_type, norbs
export add_orbital!, add_hopping!

include("tbstructs.jl")

end #module
