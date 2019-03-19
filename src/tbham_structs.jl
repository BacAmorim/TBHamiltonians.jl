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
    sorted::Vector{Bool}
end