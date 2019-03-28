## utilities

"""
    dimension(sys::BlochHamiltonian{D})  where {D}

Returns the spacial dimensional of a tight-binding Hamiltonian.
"""
function dimension(sys::BlochHamiltonian{D})  where {D}
    return D
end

"""
    norbs(sys::BlochHamiltonian)

Returns the number of orbitals within a unit cell of a BlochHamiltonian.
"""
function norbs(sys::BlochHamiltonian)
    return length(sys.orbitals)
end


## Initialize Hamiltonian
"""
    BlochHamiltonian(A, orbitals::Vector{Orbital}, hamk!::Function)

Given a vector of vector `A` = [a1, a2, ...], and  vector of `orbitals and a function `hamk!` , initializes a  BlochHamiltonian of dimension `D`, with no orbitals or hamiltonian function store
"""
function BlochHamiltonian(A, orbitals::Vector{Orbital}, hamk!::Function)

    return BlochHamiltonian(latticebasis(A), orbitals, hamk!)

end