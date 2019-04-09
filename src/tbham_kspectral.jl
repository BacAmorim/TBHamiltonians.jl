"""
    bands_full(sys::TBHamiltonian{D, T}, kpath::KPath{D}) where {D, T}

Evaluate the band structure of `sys` along `kpath`, using full diagonalization.
"""
function bands_full(sys::TBHamiltonian{D, T}, kpath::KPath{D}) where {D, T}
    
    res = pmap(k->eigen(Hermitian(hamiltoniank_full(sys, k))), kpath.kpts)
    empty = Vector{Array{1, ComplexF64}}(undef, 0)
    
    return KBandAmplitudes(
                            kpath, 
                            Vector{Array{1, ComplexF64}}(undef, 0), 
                            [res[k].values for k=1:kpath.nkpts]
                        )

end


"""
    bands_sparse(sys::TBHamiltonian{D, T}, kpath::KPath{D}; nbands, arround, ishift=0.1, kw...) where {D, T}

Evaluate the `nbands` of the band structure of `sys` closest to the energy `arround` for the `kpath`, using ARPACK.jl function `eigs`, using the shift and invert method. `ishift` is a shift in the imaginary direction to avoid poles in the real line.
"""
function bands_sparse(sys::TBHamiltonian{D, T}, kpath::KPath{D}; nbands, arround, ishift=0.1, kw...) where {D, T}
    
    res = pmap(k->eigs(hamiltoniank_sparse(sys, k); nev=nbands, which=:LM, sigma = arround + ishift*im, kw...), kpath.kpts)
    empty = Vector{Array{1, ComplexF64}}(undef, 0)
    
    return KBandAmplitudes(
                            kpath, 
                            Vector{Array{1, ComplexF64}}(undef, 0), 
                            [real.(res[k][1]) for k=1:kpath.nkpts]
                        )

end
# add dos_full and dos_sparse
# add arpes_full and arpes_sparse
# add pdos_full and pdos_sparse