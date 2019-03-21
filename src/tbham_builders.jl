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

    return TBHamiltonian(basis, Orbital[], Hopping{length(A), hoppingtype}[], [false])

end



## Add orbitals
"""
    add_orbital!(sys::TBHamiltonian, pos::Vector; l::Int=0, m::Int=0)

Given a label, position (pos), angular momentum (l, m) adds an Orbital(pos, (l, m)) to the tight-binding TBHamiltonian
"""
function add_orbital!(sys::TBHamiltonian, pos::Vector; l::Int=0, m::Int=0)


    push!(sys.orbitals, Orbital(pos, l=l, m=m))
end


"""
    add_orbital!(sys::TBHamiltonian, pos::Vector, wf::Symbol)

Given a label, position (pos), atomic wavefunction type adds an Orbital(pos, wf) to the tight-binding TBHamiltonian
"""
function add_orbital!(sys::TBHamiltonian, pos::Vector, wf::Symbol)


    push!(sys.orbitals, Orbital(pos, wf))
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

    push!(sys.hoppings, 
        Hopping(
            from, to, SVector{dimension(sys)}(sdir), convert(hopping_type(sys), val))
        )
end

