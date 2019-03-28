struct EigenStates{T}
    states::Matrix{T}
    energies::Vector{Float64}
end


struct KEigenStates{D, T}
    kmesh::KMesh{D}
    states::Vector{Matrix{T}}
    energies::Vector{Vector{Float64}}
end

struct KEigenAmplitudes{D, T, N}
    kmesh::KMesh{D}
    amplitudes::Vector{Array{N, T}}
    energies::Vector{Vector{Float64}}
end