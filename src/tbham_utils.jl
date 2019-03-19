# utility functions to extract information from TBHamiltonian
"""
    dimension(sys::TBHamiltonian{D, T})  where {D, T}

Returns the spacial dimensional of a tight-binding Hamiltonian.
"""
function dimension(sys::TBHamiltonian{D, T}) where {D, T}
    return D
end

"""
    hopping_type(sys::TBHamiltonian{D, T}) where {D, T}

Returns the type of hopping of a tight-binding Hamiltonian.
"""
function hopping_type(sys::TBHamiltonian{D, T})  where {D, T}
    return T
end

"""
    norbs(sys::TBHamiltonian)

Returns the number of orbitals within a unit cell of a tight-binding TBHamiltonian.
"""
function norbs(sys::TBHamiltonian)
    return length(sys.orbitals)
end


"""
    sort_hoppings!(sys::TBHamiltonian)

Sorts the hoppings of the TBHamiltonian (sys.hoppings), such that they are ordered first by the start orbital, then by the destination orbital and finally by the direction indices. This is necessary to work with sparce Hamiltonians.
"""
function sort_hoppings!(sys::TBHamiltonian)

    sort!(sys.hoppings, by = hop -> (hop.from, hop.to, hop.dir...))

end