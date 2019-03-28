"""
    toSVector(v::AbstractVector, len::Int)

Convert a Vector `v` to a StaticVector of length `D`. If the length of `v` is larger than `len`, the vector is truncated, otherwise the vector is padded to the right with zeros. 
"""
function toSVector(v::AbstractVector, len::Int)

    u = zeros(eltype(v), len)
    
    for i=1:min(len, length(v))
        u[i] = v[i]
    end
    
    return SVector{len}(u)
end

"""
    toSVector(v::AbstractVector)

Convert a Vector `v` to a StaticVector.
"""
toSVector(v::AbstractVector) = SVector{length(v)}(v)


"""
    cdot(a::StaticVector{d1, T1}, b::StaticVector{d2, T2}) where {d1, d2, T1, T2}

Take the inner product of two vectors of unequal size. The longer vector is truncated to the size of the smaller one. 
"""
@generated function cdot(a::StaticVector{d1, T1}, b::StaticVector{d2, T2}) where {d1, d2, T1, T2}
    s = min(d1, d2)
    
    expr = :(a[1]*b[1])
    
    for j = 2:s
        expr = :($expr + a[$j] * b[$j])
    end
    
    return :($expr)
end

"""
    cdot(a::AbstractVector, b::AbsctractVector)

Take the inner product of two vectors of unequal size. The longer vector is truncated to the size of the smaller one. 
"""
function cdot(a::AbstractVector, b::AbstractVector) 
    s = min(length(a), length(b))
    
    accu = a[1]*b[1]
    
    for j = 2:s
        accu += accu + a[j]*b[j]
    end
    
    return accu
end

