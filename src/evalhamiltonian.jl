# build hamiltonian as a function of Bloch momentum
"""
    hamiltoniank_full!(mat::AbstractMatrix, sys::TBHamiltonian{D, T1}, k::SVector{D, T2}) where {D, T1<:Number, T2<:Number}

Returns the Bloch Hamiltonian (as a full matrix) of a periodic system, `sys` as a function of the momentum `k`.
The hoppings of the Hamiltonian must be sorted
"""
function hamiltoniank_full!(mat::AbstractMatrix, sys::TBHamiltonian{D, T1}, k::SVector{D, T2}) where {D, T1<:Number, T2<:Number}

    @assert size(mat) == (D, D)

    i, j = sys.hoppings[1].to, sys.hoppings[1].from
    accu = 0.0im

    for hop in sys.hoppings


        if (i, j) == (hop.to, hop.from)
            # while from and to are the same accumulate sum_n h^{i,j}_{n} e^{ikn}
            phase = cdot(k, sys.basis*hop.dir) + cdot(k, sys.orbitals[hop.to].pos - sys.orbitals[hop.from].pos) 

            accu += hop.val*exp(-1im*phase)

        else
            # from and to have changed: store from, to and value of sum_n h^{i,j}_{n} e^{ikn}
            mat[i, j] = accu

            # reset accumulator
            phase = cdot(k, sys.basis*hop.dir) + cdot(k, sys.orbitals[hop.to].pos - sys.orbitals[hop.from].pos) 

            accu = hop.val*exp(-1im*phase)

            # update i, j
            i, j = hop.to, hop.from

        end


    end

    # write last entry of matrix
    mat[i, j] = accu


    return mat

end


"""
    hamiltoniank_full(sys::TBHamiltonian{D, T1}, k::SVector{D, T2}) where {D, T1<:Number, T2<:Number}

Returns the Bloch Hamiltonian (as a full matrix) of a periodic system, `sys` as a function of the momentum `k`.
The hoppings of the Hamiltonian must be sorted
"""
function hamiltoniank_full(sys::TBHamiltonian{D, T1}, k::SVector{D, T2}) where {D, T1<:Number, T2<:Number}

    dim = norbs(sys)
    mat = zeros(Complex{Float64}, dim, dim)

    hamiltoniank_full!(mat, sys, k)

    return mat
end


"""
    hamiltoniank_sparse(sys::TBHamiltonian{D, T1}, k::SVector{D, T2}) where {D, T1<:Number, T2<:Number}

Returns the Bloch Hamiltonian (as a sparse matrix) of a periodic system, `sys` as a function of the momentum `k`.
The hoppings of the Hamiltonian must be sorted
"""
function hamiltoniank_sparse(sys::TBHamiltonian{D, T1}, k::SVector{D, T2}) where {D, T1<:Number, T2<:Number}

    iA = Int[]
    jA = Int[]
    V = Complex{Float64}[]

    i, j = sys.hoppings[1].to, sys.hoppings[1].from
    accu = 0.0im

    for hop in sys.hoppings


        if (i, j) == (hop.to, hop.from)
            # while from and to are the same accumulate sum_n h^{i,j}_{n} e^{ikn}
            phase = cdot(k, sys.basis*hop.dir) + cdot(k, sys.orbitals[hop.to].pos - sys.orbitals[hop.from].pos) 

            accu += hop.val*exp(-1im*phase)

        else
            # from and to have changed: store from, to and value of sum_n h^{i,j}_{n} e^{ikn}
            push!(iA, i)
            push!(jA, j)
            push!(V, accu)

            # reset accumulator
            phase = cdot(k, sys.basis*hop.dir) + cdot(k, sys.orbitals[hop.to].pos - sys.orbitals[hop.from].pos) 
        
            accu = hop.val*exp(-1im*phase)

            # update i, j
            i, j = hop.to, hop.from

        end


    end

    # write last entry of matrix
    push!(iA, i)
    push!(jA, j)
    push!(V, accu)


    return sparse(iA, jA, V)

end

"""
    hamiltoniank(sys::TBHamiltonian{D, T1}, k::SVector{D, T2}; sparse=false) where {D, T1<:Number, T2<:Number}

Returns the Bloch Hamiltonian (as a sparse or full matrix) of a periodic system, `sys` as a function of the momentum `k`.
The hoppings of the Hamiltonian must be sorted.
"""
function hamiltoniank(sys::TBHamiltonian{D, T1}, k::SVector{D, T2}; sparse=false) where {D, T1<:Number, T2<:Number}
    
    if sparse
        return hamiltoniank_sparse(sys, k)
    else
        return hamiltoniank_full(sys, k)
    end
end

"""
    hamiltoniank(sys::TBHamiltonian{D, T}, k::AbstractVector; sparse=false) where {D, T<:Number}

Returns the Bloch Hamiltonian (as a sparse or full matrix) of a periodic system, `sys` as a function of the momentum `k`.
The hoppings of the Hamiltonian must be sorted. 
It is checked if the hoppings of the hamiltonian are sorted and if the length of the `k` vectors matches the dimensionality of the system.
"""
function hamiltoniank(sys::TBHamiltonian{D, T}, k::AbstractVector; sparse=false) where {D, T<:Number}
    @assert length(k) == D "length of k vector much be the same as the dimensionality of the system"
    @assert sys.sorted[1] "hoppings of the Hamiltonian must be sorted"
    
    return hamiltoniank(sys, k, sparse=sparse)
end


"""
    hamiltoniank!(mat::AbstractMatrix, sys::TBHamiltonian{D, T}, k::AbstractVector; sparse=false) where {D, T<:Number}

In-place version of `hamiltoniank`. Only works for full matrices.
"""
function hamiltoniank!(mat::AbstractMatrix, sys::TBHamiltonian{D, T}, k::AbstractVector) where {D, T<:Number}
    @assert length(k) == D "length of k vector much be the same as the dimensionality of the system"
    @assert sys.sorted[1] "hoppings of the Hamiltonian must be sorted"

    hamiltoniank_full!(mat, sys, k)

end