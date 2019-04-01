function R2d(theta) = SMatrix{2, 2}([cos(theta) -sin(theta); sin(theta) cos(theta)])

function tblg_commensurate_angle(m, r)

    if gcd(m, r) != 1
        println("m and r are not coprime")
        nothing
    end

    return acos((3*m^2 + 3*m*r + r^2/2)/(3*m^2 + 3*m*r + r^2))
end


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

function tblg_commensurate_number_atoms(m, r)

    if gcd(r, 3) == 1
        return 4*(3*m^2 + 3*m*r + r^2)
    elseif gcd(r, 3) == 3
        return 4*(m^2 + m*r + div(r^2, 3))
    end

end

function tblg_commensurate_carbon_positions(m, r, a0)

    @assert gcd(m, r) ==1 "m and r must be coprimes"

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
    τAbot = R2d(pi/2-α)*[0.0, 0.0]
    τBbot = R2d(pi/2-α)*[0.0, -a0/sqrt(3)]

    # bottom layer basis vectors
    a1top = R2d(θ + pi/2 - α)*a1
    a2top = R2d(θ + pi/2 - α)*a2
    τAtop = R2d(θ + pi/2-α)*SVector(0.0, 0.0)
    τBtop = R2d(θ + pi/2-α)*SVector(0.0, a0/sqrt(3))

    N = ceil(Int, norm(A1)/a0)
    listAbot = SVector{Float64}[]
    listBbot = SVector{Float64}[]
    listAtop = SVector{Float64}[]
    listBtop = SVector{Float64}[]

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



function tblg_commensurate_build(m, r, a0, d, t, Vinter::Function, hoprange)

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