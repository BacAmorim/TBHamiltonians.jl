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