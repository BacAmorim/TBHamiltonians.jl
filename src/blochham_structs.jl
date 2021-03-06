# added structure for Bloch Hamiltonian. Instead of storing  a list of vectors, an in-place function is stored.

"""
    BlochHamiltonian{D, F<:Function}

Stores a tight-binding system in Bloch space.  

### Fields
The tight-binsing system consists of:
- 'latbasis::SMatrix{D, D, Float64}': the latbasis for a `D`-dimensional Bravais lattice. The basis vectors are stored as columns of the matrix basis, i. e. a1 = basis[:, 1]
- 'orbitals::Vector{Orbital}': a list of orbitals within the unit cell
- 'hk!::F': The Bloch Hamiltonian: a function that for a given `k` acts on a matrix `mat` and stores there H(k)
"""
struct BlochHamiltonian{D, F}
    basis::SMatrix{D, D, Float64}
    orbitals::Vector{Orbital}
    hk!::F
end

"""
    BlochHamiltonian_noF{D}

Stores a tight-binding system in Bloch space.  

### Fields
The tight-binsing system consists of:
- 'latbasis::SMatrix{D, D, Float64}': the latbasis for a `D`-dimensional Bravais lattice. The basis vectors are stored as columns of the matrix basis, i. e. a1 = basis[:, 1]
- 'orbitals::Vector{Orbital}': a list of orbitals within the unit cell
- 'hk!::F': The Bloch Hamiltonian: a function that for a given `k` acts on a matrix `mat` and stores there H(k)
"""
struct BlochHamiltonian_noF{D}
    basis::SMatrix{D, D, Float64}
    orbitals::Vector{Orbital}
    hk!::Function
end