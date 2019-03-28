"""
    reciprocal(A::AbstractMatrix)

Given a basis `A = [a1 a2 ...]` (with the basis vectors forming the matrix columns), returns the reciprocal basis `B = [b1 b2 ...]` (with the reciprocal basis vectors forming the matrix columns) such that:
dot(aᵢ, bⱼ) = 2πδᵢⱼ
"""
reciprocal(A::AbstractMatrix) = 2π*transpose(inv(A))