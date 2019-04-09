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

struct KBandStates{D, T}
    kpath::KPath{D}
    states::Vector{Matrix{T}}
    energies::Vector{Vector{Float64}}
end

struct  KBandAmplitudes{D, T, N}
    kpath::KPath{D}
    amplitudes::Vector{Array{N, T}}
    energies::Vector{Vector{Float64}}
end