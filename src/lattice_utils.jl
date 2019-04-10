"""
    latticebasis(A::AbstractVector...)

Given a tuple of basis vectors `A = (a1,  a2, ...)` returns a SMatrix with the vector basis as columns.
"""
function latticebasis(A::AbstractVector...)
    dim = length(A)
    for ai in A
        @assert length(ai) == dim "number of basis vectors, must coincide with their dimension"
    end
    
    return SMatrix{dim, dim}(hcat(A...))
end

"""
    reciprocalbasis(A::AbstractMatrix)

Given a basis `A = [a1 a2 ...]` (with the basis vectors forming the matrix columns), returns the reciprocal basis `B = [b1 b2 ...]` (with the reciprocal basis vectors forming the matrix columns) such that:
dot(aᵢ, bⱼ) = 2πδᵢⱼ
"""
reciprocalbasis(A::AbstractMatrix) = 2π*transpose(inv(A))


"""
    spanlattice(B::SMatrix, cutoff::Float64)

Span the lattice points generated by the lattice basis vectors B[:, i], with length shorter than `cutoff`.
"""
function spanlattice(B::SMatrix, cutoff::Float64)
    
    Ns = [ceil(Int, cutoff/norm(B[:, i])) for i=1:size(B, 2)]
    dim = size(B, 1)

    itr = Iterators.product([-Ni:Ni for Ni in Ns]...)
    
    list=SVector{dim, Float64}[]
    
    for I in itr
        
        pt = B*SVector(I...)
        if norm(pt) <= cutoff
            push!(list, pt)
        end
    end
    
    return list
end


"""
    addlattice(lats::Vector{Vector{StaticVectors}}, cutoff::Float64)

Given n spanned lattices [lats[1], ..., lats[N]], form tuples of vector of the form (R₁, ..., Rₙ) with Rᵢ belonging to lats[i],  such that norm(R₁ + ... + Rₙ) < cutoff.
"""
function addlattice(lats::Vector{Vector{SVector{D, T}}}, cutoff::Float64) where {D, T}
    
    itr = Iterators.product(lats...)
    Nlat = length(lats)
    list = NTuple{Nlat, SVector{D, T}}[]
    
    for G in itr
        if norm(sum(G)) < cutoff
            push!(list, G)
        end
    end
    
    return list
end