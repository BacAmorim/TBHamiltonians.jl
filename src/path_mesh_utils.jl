"""
    path(vertices, npts)

Creates a path in k-space. The path is formed by lines connecting the points provided in `vertices` as an array of arrays.
The return path has approx `npts` points.
The path is returned as a `KPath` object.
"""
function path(vertices, npts)

    for pt in vertices
        @assert length(pt)==length(vertices[1]) "all points must be of the same size"
    end

    v = [toSVector(u) for u in vertices]

    lengths = [norm(v[i+1] - v[i]) for i=1:length(v)-1]
    plengths = [sum(lengths[1:i]) for i=1:length(lengths)]

    totallength = plengths[end]

    dx = totallength/npts

    pNpts = floor.(lengths/dx)

    kpts = [v[1] + (v[2]-v[1])*(n-1)/pNpts[1] for n=1:pNpts[1]]
    arc = [norm((v[2]-v[1])*(n-1)/pNpts[1]) for n=1:pNpts[1]]

    for p=2:length(vertices)-1

        append!(kpts, [v[p] + (v[p+1]-v[p])*(n-1)/pNpts[p] for n=1:pNpts[p]])
        append!(arc, 
                [
                arc[end] + norm((v[p]-v[p-1])/pNpts[p-1]) + norm((v[p+1]-v[p])*(n-1)/pNpts[p]) for n=1:pNpts[p]]
            )
    end

    append!(kpts, [v[end]])
    append!(arc, [arc[end] + norm((v[end]-v[end-1])/pNpts[end])])
    prepend!(plengths, [0.0])


    return KPath(kpts, arc, plengths, length(kpts))

end


"""
    mpmesh(B::AbstractMatrix, Ns::NTuple{D, Int}) where {D}

Creates a Ns[1]*Ns[2]*... Monkhorst-Pack mesh from the basis vectors B[:, 1], B[:, 2], ...
The mesh is returned as a `KMesh` object.
"""
function mpmesh(B::AbstractMatrix, Ns::NTuple{D, Int}) where {D}

    @assert size(B, 1) == size(B, 2) "number of basis vectors must be the same as a the dimension of the vectors"
    @assert D == size(B, 2) "dimensions of mesh must coincide with the number of basis vectors"
    
    itr = Iterators.product([1:Ni for Ni in Ns ]...)
    
    sB = SMatrix{size(B)...}(B)
    
    list = SVector[]
    for I in itr
        push!(list, sB*SVector(I...))
    end
    
    len = length(list)

    return KMesh(list, [1/len for k in list], len)
end