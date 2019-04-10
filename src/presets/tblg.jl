"""
    R2d(theta)

Two-dimensional rotation matrix by an angle `theta`. 
"""
function R2d(theta) 
    SMatrix{2, 2}([cos(theta) -sin(theta); sin(theta) cos(theta)])
end

"""
    tblg_commensurate_angle(m, r)

Returns the commensurate angle for a commensurate twisted bilayer graphene structure labelled by co-prime integers `m` and `r`, according to the notation of PRB 86, 155449 (2012) [arXiv:1202.1088].
"""
function tblg_commensurate_angle(m, r)

    if gcd(m, r) != 1
        println("m and r are not coprime")
        nothing
    end

    return acos((3*m^2 + 3*m*r + r^2/2)/(3*m^2 + 3*m*r + r^2))
end

"""
    tblg_commensurate_basis(m, r, a0)

Returns the Bravais basis for a commensurate (m, r) twisted bilayer graphene structure. `a0` is the lattice parameter of an isolated graphene layer.
"""
function tblg_commensurate_basis(m, r, a0)

    a1 = a0*SVector(0.5, sqrt(3)/2)
    a2 = a0*SVector(-0.5, sqrt(3)/2)

    if gcd(r, 3) == 1
        t1 = (m)*a1 + (m+r)*a2
        t2 = -(m+r)*a1 + (2*m+r)*a2

    elseif gcd(r, 3) == 3
        t1 = (m+div(r,3))*a1 + div(r, 3)*a2
        t2 = -div(r, 3)*a1 + (m+2*div(r, 3))*a2
    end

    return t1, t2

end

"""
    tblg_commensurate_number_atoms(m, r)

Number of total carbon atoms in a unit cell of a commensurate (m, r) twisted bilayer graphene structure.
"""
function tblg_commensurate_number_atoms(m, r)

    if gcd(r, 3) == 1
        return 4*(3*m^2 + 3*m*r + r^2)
    elseif gcd(r, 3) == 3
        return 4*(m^2 + m*r + div(r^2, 3))
    end

end

"""
    tblg_commensurate_carbon_positions(m, r, a0)

Determines the position of the carbon atoms inside a unit cell of a commensurate (m, r) twisted bilayer graphene structure. `a0` is the lattice parameter of a single layer graphene. The functions returns as a tuple:
θ, Natls, WZ, [A1, A2], [a1bot, a2bot], [a1top, a2top], [listAbot, listBbot, listAtop, listBtop]
where:
- θ is the commensurate twist angle,
- Natls is number of carbon atoms per layer, per sublattice in the unit cell, 
- WZ is a vector of positions, deliminting the Wigner-Seitz unit cell of the tBLG structure,
- [A1, A2] is the Bravais basis of the commensurate structure 
- [a1bot, a2bot] is the Bravais basis of the isolated bottom layer
- [a1top, a2top] is the Bravais basis of the isolated top layer
- [listAbot, listBbot, listAtop, listBtop] constains the position of the atoms in each layer and sublattice. 
"""
function tblg_commensurate_carbon_positions(m::Int, r::Int, a0)

    # @assert gcd(m, r) == 1 "m and r must be coprimes"

    θ = tblg_commensurate_angle(m, r)
    Natls =div(tblg_commensurate_number_atoms(m, r), 4)

    a1 = a0*SVector(0.5, sqrt(3)/2)
    a2 = a0*SVector(-0.5, sqrt(3)/2)

    A1, A2 = tblg_commensurate_basis(m, r, a0)

    α = atan(A1[2]+A2[2], A1[1]+A2[1])

    #rotate basis vectors
    A1 = R2d(pi/2-α)*A1
    A2 = R2d(pi/2-α)*A2

    # WZ unit cell
    C1 = (A1 + A2)/3
    WZ = [R2d(2*pi*n/6)*C1 for n=0:5]

    # bottom layer basis vectors
    a1bot = R2d(pi/2-α)*a1
    a2bot = R2d(pi/2-α)*a2
    τAbot = R2d(pi/2-α)*SVector(0.0, 0.0)
    τBbot = R2d(pi/2-α)*SVector(0.0, -a0/sqrt(3))

    # top layer basis vectors
    a1top = R2d(θ + pi/2 - α)*a1
    a2top = R2d(θ + pi/2 - α)*a2
    τAtop = R2d(θ + pi/2-α)*SVector(0.0, 0.0)
    τBtop = R2d(θ + pi/2-α)*SVector(0.0, a0/sqrt(3))

    N = ceil(Int, norm(A1)/a0)
    listAbot = SVector{2, Float64}[]
    listBbot = SVector{2, Float64}[]
    listAtop = SVector{2, Float64}[]
    listBtop = SVector{2, Float64}[]

    for n=-N:N, m=-N:N
        pbot = n*a1bot + m*a2bot
        if isinsideconvex(pbot, WZ)
            push!(listAbot, pbot + τAbot)
            push!(listBbot, pbot + τBbot)
        end

        ptop = n*a1top + m*a2top
        if isinsideconvex(ptop, WZ)
            push!(listAtop, ptop + τAtop)
            push!(listBtop, ptop + τBtop)
        end

    end

    # test we captured all the points
    if length(listAbot) != Natls
        println("Warning: Not all A atoms from bottom layer were found")
    end
    if length(listBbot) != Natls
        println("Warning: Not all B atoms from bottom layer were found")
    end
    if length(listAtop) != Natls
        println("Warning: Not all A atoms from top layer were found")
    end
    if length(listBtop) != Natls
        println("Warning: Not all B atoms from top layer were found")
    end

    return θ, Natls, WZ, [A1, A2], [a1bot, a2bot], [a1top, a2top], [listAbot, listBbot, listAtop, listBtop]

end


"""
    Vinter_tblg_sk(x; Vppp0 = -2.7, Vpps0 = 0.48, a0 = 2.46, d = 3.35, r0 = 0.45264)

Slater-Koster parametrization of the interlayer hopping in twisted bilayer graphene as a function of the in-plane separation `x` (a vector).
"""
function Vinter_tblg_sk(x; Vppp0 = -2.7, Vpps0 = 0.48, a0 = 2.46, d = 3.35, r0 = 0.45264)

    r = norm(x)

    return Vpps0*d^2/(r^2+d^2)*exp(-(sqrt(r^2+d^2)-d)/r0) + Vppp0*r^2/(r^2+d^2)*exp(-(sqrt(r^2+d^2)-a0/sqrt(3))/r0)

end


"""
    tblg_commensurate_build(m, r, a0, d, t, hoprange, Vinter::Function=Vinter_tblg_sk)

Builds a TBHamiltonian system for a commensurate (m, r) twisted bilayer graphene structure. The remaining arguments are:
- `a0` is the lattice parameter of single layer graphene;
- `d` is the interlayer seperation;
- `t` is the nearest-neighbour in-plane hopping;
- `hoprange` is a distance cutoff for the interlayer hoppings;
- `Vinter` is a function parametrizing the interlayer hopping as a function of the in-plane seperation.

It resturns: tblg, theta, Napl, WZ, abot, atop; where:
- `tblg`, is the TBHamiltonian struct;
- `theta`, is the twist angle (in radians);
- `Napl`, is the number of atoms per layer and per sublattice in the unit cell;
- `WZ`, are the vertices of the Wigner-Seitz unit cell;
- `abot`, is the Bravais basis of the botton layer;
- `atop`, is the Bravais basis of the top layer.
"""
function tblg_commensurate_build(m, r, a0, d, t, hoprange, Vinter::Function=Vinter_tblg_sk)

    theta, Napl, WZ, A, abot, atop, carbon = tblg_commensurate_carbon_positions(m, r, a0)

    tblg = TBHamiltonian(A[1], A[2])

    # Add orbtials
    ## A top atoms
    for atom in carbon[3]
        add_orbital!(tblg, SVector(atom..., d), :pz)
    end

    ## B top atoms
    for atom in carbon[4]
        add_orbital!(tblg, SVector(atom..., d), :pz)
    end

    ## A bot atoms
    for atom in carbon[1]
        add_orbital!(tblg, SVector(atom..., 0.0), :pz)
    end

    ## B bot atoms
    for atom in carbon[2]
        add_orbital!(tblg, SVector(atom..., 0.0), :pz)
    end


    # Add intralayer hoppings
    ## nearest-neighbours unit cells
    NN = [(0, 0), (1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]

    ## create trees for NN atom search
    topR = zeros(2, 2*Napl)
    for i=1:Napl
        topR[:, i] = carbon[3][i]
        topR[:, i + Napl] = carbon[4][i]
    end
    toptree = BallTree(topR)

    botR = zeros(2, 2*Napl)
    for i=1:Napl
        botR[:, i] = carbon[1][i]
        botR[:, i + Napl] = carbon[2][i]
    end
    bottree = BallTree(botR)

    ## add intralayer hopping in top layer
    for (n, m) in NN
        for i=1:2*Napl
            idxs = inrange(toptree, n*A[1] + m*A[2] + topR[:, i], 1.1*a0/sqrt(3), true)
            for j in idxs
                if norm(n*A[1] + m*A[2] + topR[:, i] - topR[:, j])>0.1*a0/sqrt(3)
                    add_hopping!(tblg, j, i, [n, m], -t)
                end
            end
        end
    end

    ## add intralayer hopping in bottom layer
    for (n, m) in NN
        for i=1:2*Napl
            idxs = inrange(bottree, n*A[1] + m*A[2] + botR[:, i], 1.1*a0/sqrt(3), true)
            for j in idxs
                if norm(n*A[1] + m*A[2] + botR[:, i] - botR[:, j])>0.1*a0/sqrt(3)
                    add_hopping!(tblg, 2*Napl + j, 2*Napl + i, [n, m], -t)
                end
            end
        end
    end

    ## add interlayer hoppings
    for (n, m) in NN
        ### from bottom to top
        for i=1:2*Napl
            idxs = inrange(bottree, n*A[1] + m*A[2] + topR[:, i], hoprange, true)
            for j in idxs
                add_hopping!(tblg, 2*Napl + j, i, [n, m], Vinter(n*A[1] + m*A[2] + topR[:, i] - botR[:, j]))
            end
        end

        ### from top to bottom
        for i=1:2*Napl
            idxs = inrange(toptree, n*A[1] + m*A[2] + botR[:, i], hoprange, true)
            for j in idxs
                add_hopping!(tblg, j, 2*Napl + i, [n, m], Vinter(n*A[1] + m*A[2] + botR[:, i] - topR[:, j]))
            end
        end
    end

    sort_hoppings!(tblg)

    return tblg, theta, Napl, WZ, abot, atop

end