## utilities

"""
    dimension(sys::BlochHamiltonian{D, F})  where {D, T}

Returns the spacial dimensional of a tight-binding Hamiltonian.
"""
function dimension(sys::BlochHamiltonian{D, F})  where {D, T}
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
    BlochHamiltonian(A, orbitals::Vector{Orbital}, F::Function)

Given a vector of vector `A` = [a1, a2, ...], and  vector of `orbitals and a function `hamk!` , initializes a  BlochHamiltonian of dimension `D`, with no orbitals or hamiltonian function store
"""
function BlochHamiltonian(A, orbitals::Vector{Orbital}, F::Function)

    for a in A
        @assert length(A)==length(a) "number of basis vectors and their dimensions must coincide"
    end
    dim = length(A)
    
    basis = SMatrix{dim, dim, Float64}(hcat(A...))

    return BlochHamiltonian(basis, orbitals, hamk!)

end

# eval hamiltonian
function hamiltoniank!(mat::AbstractMatrix, sys::BlochHamiltonian, k::AbstractVector)
    @assert length(k) == D "length of k vector much be the same as the dimensionality of the system"

    sys.hk!(mat, SVector{D}(k))

end

function hamiltoniank(sys::BlochHamiltonian, k::AbstractVector)
    @assert length(k) == D "length of k vector much be the same as the dimensionality of the system"

    dim = norbs(sys)
    mat = zeros(Complex{Float64}, dim, dim)

    sys.hk!(mat, SVector{D}(k))

end