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


# add bands_sparse
# add dos_full and dos_sparse
# add arpes_full and arpes_sparse
# add pdos_full and pdos_sparse