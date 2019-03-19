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

    return TBHamiltonian(basis, Orbital[], Hopping{length(A), hoppingtype}[], false)

end


## Add orbitals
"""
    add_orbital!(sys::TBHamiltonian, pos::Vector; l::Int=0, m::Int=0)

Given a label, position (pos), angular momentum (l, m) adds an Orbital(pos, (l, m)) to the tight-binding TBHamiltonian
"""
function add_orbital!(sys::TBHamiltonian, pos::Vector; l::Int=0, m::Int=0)

    pos3 = zeros(Float64, 3)  # conver the position vector to a 3D vector
    for i=1:min(3, length(pos))
        pos3[i] = convert(Float64, pos[i])
    end

    push!(sys.orbitals, Orbital(SVector{3}(pos3), (l, m)))
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

    push!(sys.orbitals, Orbital(SVector{3}(pos3), replacement_wf[wf]))
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

