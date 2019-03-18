# Basic structures

"""
    Orbital

Stores an orbital.

### Fields
- 'pos::SVector{3,Float64}': position of the orbital within the unit cell
- 'lm::NTuple{2, Int}': orbital quantum numbers: (total angular momentum, azimuthal angular momentum)
"""
struct Orbital
    pos::SVector{3,Float64}
    lm::NTuple{2, Int}
end

"""
    Hopping{D, T<:Number}

Stores a hopping with amplitude of type `T` in a `D`-dimensional lattice.

### Fields
- 'from::Int': start orbital
- 'to::Int': end orbital
- 'dir::SVector{D, Int}': direction of the hopping (in units of Bravais basis vectors. 
- 'val::T': value of hopping
"""
struct Hopping{D, T<:Number}
    from::Int
    to::Int
    dir::SVector{D, Int}
    val::T
end

"""
    TBHamiltonian{D, T<:Number}

Stores a tight-binding system define in real space for a `D`-dimensional lattice (embedded in 3D) with hoppings of type `T`

### Fields
The tight-binsing system consists of:
- 'basis::SMatrix{D, D, Float64}': the basis for a `D`-dimensional Bravais lattice. The basis vectors are stored as columns of the matrix basis, i. e. a1 = basis[:, 1]
- 'orbitals::Vector{Orbital}': a list of orbitals within the unit cell
- 'hoppings::Vector{Hopping{D, T}': a list of hoppings (which includ intra- an intercell hoppings)
- 'sorted::Bool': whether the list of hoppings is been sorted by (from, to). Use of sparce matrices requires sorting.
"""
struct TBHamiltonian{D, T<:Number}
    basis::SMatrix{D, D, Float64}
    orbitals::Vector{Orbital}
    hoppings::Vector{Hopping{D, T}}
    sorted::Bool
end


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

# functions for construction of Hamiltonians


## Initialize Hamiltonian
"""
    TBHamiltonian(A::Vector...; hoppingtype=Float64)

Given a tuple of vectors a1, a2, ..., initializes a tight-binding system of dimension `D`, and hoppings of type `hoppingtype`, without any orbitals or hoppings added. 
"""
function TBHamiltonian(A::Vector...; hoppingtype=Float64)

    for a in A
        @assert length(A)==length(a) "number of basis vectors and their dimensions must coincide"
    end
    dim = length(A)
    
    basis = SMatrix{dim, dim, Float64}(hcat(A...))

    return TBHamiltonian(basis, Orbital[], Hopping{length(A), hoptype}[], false)

end


## Add orbitals
"""
    add_orbital!(sys::TBHamiltonian, pos::Vector, l::Int=0, m::Int=0)

Given a label, position (pos), angular momentum (l, m) adds an Orbital(pos, (l, m)) to the tight-binding TBHamiltonian
"""
function add_orbital!(sys::TBHamiltonian, pos::Vector; l::Int=0, m::Int=0)

    pos3 = zeros(Float64, 3)  # conver the position vector to a 3D vector
    for i=1:min(3, length(pos))
        pos3[i] = convert(Float64, pos[i])
    end

    add_orbital!(sys, Orbital(SVector{3}(pos3), (l, m)))
end

replacement_wf = Dict(
    :s => (0, 0),
    :px => (1, 1),
    :py => (1, -1),
    :pz => (1, 0),
    :dxy => (2, -2),
    :dx2y2 => (2, 2),
    :dxz => (2, 1),
    :dyz => (2, -1),
    :dz2 => (2, 0)
)

"""
    add_orbital!(sys::TBHamiltonian, pos::Vector, wf::Symbol)

Given a label, position (pos), atomic wavefunction type adds an Orbital(pos, wf) to the tight-binding TBHamiltonian
"""
function add_orbital!(sys::TBHamiltonian, pos::Vector, wf::Symbol)

    @assert wf in keys(replacement_wf) "Shorthand for wavefunction type unknown. Available shorthands: `:s`, `:px`, `:py`,`:pz`, `:dxy`, `:dx2y2`, `:dxz`, `:dyz`, `:dz2`."

    pos3 = zeros(Float64, 3)  # conver the position vector to a 3D vector
    for i=1:min(3, length(pos))
        pos3[i] = convert(Float64, pos[i])
    end

    add_orbital!(sys, Orbital(SVector{3}(pos3), replacement_wf[wf]))
end


## Add hopping
"""
    add_hopping!(sys::TBHamiltonian, from::Int, to::Int, dir::Vector{Int}, val)

Given a initial orbital index (from), final orbital (to), direction (dir) and a value (val), adds a Hopping(dir, from, to, val) to the tight-binding system.
"""
function add_hopping!(sys::TBHamiltonian, from::Int, to::Int, dir::Vector{Int}, val)

    sdir = zeros(Int, dimension(sys)) # conver the direction vector to a 3D vector
    for i=1:min(dimension(sys), length(dir))
        sdir[i] = dir[i]
    end

    push!(sys.hoppings, Hopping(SVector{dimension(sys)}(sdir), 
            from, to, 
            convert(hopping_type(sys), val))
        )
end