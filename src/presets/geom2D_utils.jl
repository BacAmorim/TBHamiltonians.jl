"""
    triangle_area(v1, v2, v3)

Area of a triangle defined by the vertices `v1`, `v2` and `v3`.
"""
function triangle_area(v1, v2, v3)
    u1 = v2 - v1
    u2 = v3 - v1

    return abs(u1[1]*u2[2] - u1[2]*u2[1])/2

end

"""
    isinsideconvex(p, vertices; tol=sqrt(eps()))

Determines if a point `p` is inside a convex region defined by the vertices in the array `vertices`.
"""
function isinsideconvex(p, vertices; tol=sqrt(eps()))
    n = length(vertices)

    center = [0.0, 0.0]
    for i=1:n
        center += vertices[i]
    end
    center = center/n

    convex_area = 0.0
    for i=1:n
        convex_area += triangle_area(center, vertices[i], vertices[mod1(i+1,n)])
    end

    test_area = 0.0
    for i=1:n
        test_area += triangle_area(p, vertices[i], vertices[mod1(i+1,n)])
    end

    return isapprox(test_area, convex_area, rtol = tol)
end
