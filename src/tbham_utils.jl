# function to sort hoppings
"""
    sort_hoppings!(sys::TBHamiltonian)

Sorts the hoppings of the TBHamiltonian (sys.hoppings), such that they are ordered first by the start orbital, then by the destination orbital and finally by the direction indices. This is necessary to work with sparce Hamiltonians.
"""
function sort_hoppings!(sys::TBHamiltonian)

    sort!(sys.hoppings, by = hop -> (hop.from, hop.to, hop.dir...))

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
    hopping_ranges(sys::TBHamiltonian)

Returns a vector with the range of hoppings in the TBHamiltonian, i.e., for a list of hopping in 3D along the directions n*a1 + m*a2+ k*a3, returns maximum(n), maximum(m), maximum(k)
"""
function hopping_ranges(sys::TBHamiltonian)

    dim = dimension(sys)
    ranges = [[0, 0] for n=1:dim]

    for hop in sys.hoppings
        for i=1:dim
            if ranges[i][1] > hop.dir[i] #min of direction i
                ranges[i][1] = hop.dir[i]
            end
            if ranges[i][2] < hop.dir[i] # max of direction i
                ranges[i][2] = hop.dir[i]
            end
        end
    end

    return [ranges[i][1]:ranges[i][2] for i=1:dim]
end


"""
    hopping_array(sys::TBHamiltonian)

Returns the hoppings t^{to, from}_{n1, n2, ...} store in a array t[n1, n2, .., to, from]. n1, ... run from -range1, ..., range1
"""
function hopping_array(sys::TBHamiltonian)
    dim = dimension(sys)

    ranges = hopping_ranges(sys)
    Norbs = norbs(sys)

    t = OffsetArray{Float64, dim + 2}(ranges..., 1:Norbs, 1:Norbs)
    for i in eachindex(t)
        t[i] = zero(hopping_type(sys))
    end

    for hop in sys.hoppings
        t[hop.dir..., hop.to, hop.from] = hop.val
    end

    return t # array t^{t,j}
end