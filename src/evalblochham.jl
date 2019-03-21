# eval hamiltonian
function hamiltoniank!(mat::AbstractMatrix, sys::BlochHamiltonian, k::AbstractVector)
    @assert length(k) == D "length of k vector much be the same as the dimensionality of the system"

    sys.hk!(mat, SVector{D}(k))

end

function hamiltoniank(sys::BlochHamiltonian, k::AbstractVector)
    @assert length(k) == D "length of k vector much be the same as the dimensionality of the system"

    dim = norbs()
    mat = zeros(Complex{Float64}, )

    sys.hk!(mat, SVector{D}(k))

end