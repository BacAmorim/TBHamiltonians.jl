"""
    KPath{D}

Stores an path in k-space.

#### Fields
- kpts::Vector{SVector{D, Float64}}
- arc::Vector{Float64}
- special::Vector{Float64}
- nkpts::Int

"""
struct KPath{D}
    kpts::Vector{SVector{D, Float64}}
    arc::Vector{Float64}
    special::Vector{Float64}
    nkpts::Int
end

"""
    KMesh{D}

Stores an path in k-space.

#### Fields
- kpts::Vector{SVector{D, Float64}}
- ws::Vector{Float64}
- nkpts::Int
"""
struct KMesh{D}
    kpts::Vector{SVector{D, Float64}}
    ws::Vector{Float64}
    nkpts::Int
end