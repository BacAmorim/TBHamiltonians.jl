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
    Orbital(pos::Vector; l::Int=0, m::Int=0)

Initialized an orbital centered at `pos` of type (l, m)"
"""
function Orbital(pos::Vector; l::Int=0, m::Int=0)
    
    pos3 = zeros(Float64, 3)  # conver the position vector to a 3D vector
    for i=1:min(3, length(pos))
        pos3[i] = convert(Float64, pos[i])
    end

    return Orbital(SVector{3}(pos3), (l, m))
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
    Orbital(pos::Vector, wf::Symbol)

Initialized an orbital centered at `pos` of type `wf`. 
Known types: `:s`, `:px`, `:py`, `:pz`, `:dxy`, `:dx2y2`, `:dxz`, `:dyz`, `:dz2`"
"""
function Orbital(pos::Vector, wf::Symbol)
    
    @assert wf in keys(replacement_wf) "Shorthand for wavefunction of known type. Available shorthands: `:s`, `:px`, `:py`,`:pz`, `:dxy`, `:dx2y2`, `:dxz`, `:dyz`, `:dz2`."

    pos3 = zeros(Float64, 3)  # conver the position vector to a 3D vector
    for i=1:min(3, length(pos))
        pos3[i] = convert(Float64, pos[i])
    end

    return Orbital(pos, replacement_wf[:wf])
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