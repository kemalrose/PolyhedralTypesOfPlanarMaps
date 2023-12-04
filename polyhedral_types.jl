


# Authors: Kemal Rose, Boulos es Hilany
# Date: November 14, 2023
# Short description: This code accompanies the paper 
# "Lower bounds on topological types of maps" by Kemal Rose and Boulos es Hilany
import Pkg;
Pkg.add("Polymake")
Pkg.update()
using Oscar, LinearAlgebra, Random, Polymake


# Takes a 2xn matrix and checks wether the columns represent vertices of a polygon in clockwise order, starting with the lexicographically smallest one.
# -------------  Input:
# pols              an integer 2xn matrix
# -------------  Output:
# new_pols          a boolean

function isCorrect(pol)
    res = []
    nverts = size(pol,2)
    if nverts == 2
        return pol[:,1] < pol[:,2] 
    else 
        for i in 1:nverts
            s = pol[:,(i-1)%nverts+1]
            p = pol[:,((i-1)+1)%nverts+1]
            q = pol[:,((i-1)+2)%nverts+1]
            push!(res,isAdmissible(s,p,q))
        end
    end
    all(res)
end

# Lists all integer polygons contained in d times the 2-simplex
# -------------  Input:
# d              a natural number.
# -------------  Output:
# polytopeMatrix   a 2xn matrix. The columns represent vertices in clockwise order, starting with the lexicographically smallest one.

function list_all_polygons(d)
    S = [[i,j] for i in 0:d for j in 0:d]
    filter!(s-> sum(s) ≤ d && sum(s) > 0, S)
    S = sort(S)
    polytopeMatrix = []
    newPolytopes = []
    for s1 in S
        for s2 in S
            if s1 < s2
                push!(newPolytopes, hcat(s1, s2))
            end
        end
    end
    while length(newPolytopes) > 0
        polytopeMatrix = [polytopeMatrix; newPolytopes]
        newPolytopes = extendPolytopes(newPolytopes, S)
    end
    polytopeMatrix
end

# Given a matrix P representing the vertices a polygon in clockwise order, and a point v in the plane, checks wether (P∣v) also represents a polygon.
# -------------  Input:
# P              an integer 2xn matrix.
# -------------  Output:
# bool               a boolean.

function isAdmissible(P, v)
    s = P[:,1]
    l = P[:,2]
    p = P[:,end-1]
    q = P[:,end]
    if size(P,2) == 2
        return isAdmissible(s, q, v)
    else
        return isAdmissible(p, q, v) && isAdmissible(s, q, v) && isAdmissible(v, s, l) && s < v
    end
end


# Given three points p,q,v in the plane this function checks wether v lies on the right of the path p->q.
# -------------  Input:
# p,q,v              a triple of vectors with two entries each.
# -------------  Output:
# bool               a boolean.

function isAdmissible(p, q, v)
    qn = q - p
    vn = v - p
    bool = qn[2] * vn[1] - qn[1] * vn[2] > 0
    bool
end


# Given a set of matrices representing polygons and a set S of points, computes all polygons that can be obtained by adding one element of S.
# -------------  Input:
# pols              a list of integer matrices
# -------------  Output:
# new_pols          a list of integer matrices

function extendPolytopes(pols, S)
    new_pols = []
    for pol in pols
        for s in S
            if isAdmissible(pol, s)
                push!(new_pols, hcat(pol,s))
            end
        end
    end
    new_pols
end


# Checks wether a polygon is "conic".
# -------------  Input:
# verts             2xn matrix
# -------------  Output:
# bool              boolean

function is_conic(verts)
    pol = convex_hull(transpose(verts))
    lpts = Vector{Vector{Int64}}(lattice_points(pol))
    M = zeros(Int64,5,size(lpts,1))
    for j in 1:size(lpts,1)
        s1, s2 = lpts[j]
        M[1,j] = s1
        M[2,j] = s2
        M[3,j] = s1^2
        M[4,j] = s1*s2
        M[5,j] = s2^2
    end
    bool = rank(M) >= 5
    bool
end

# Lists all conic integer polygons contained in d times the 2-simplex 
# -------------  Input:
# d              a natural number.
# -------------  Output:
# polytopeMatrix   a 2xn matrix. The columns represent vertices in clockwise order, starting with the lexicographically smallest one.

function list_all_conic_pols(deg::Int)
    all_pols = list_all_polygons(deg)
    conic_pols = filter(verts->is_conic(verts), all_pols)
    conic_pols
end




# Given to polygons, computes Delta. 
# Delta is the newton polytope of the map (f1, f2): C^2 -> C^2 where f1, f2 are generic polynomials with Newton polytopes A1, A2 respectively. 
# -------------  Input:
# (A1, A2)             a pair of Polyhedra
# -------------  Output:
# Delta                a Polyhedron

function get_delta(A1, A2)
    verts1, verts2 = A1, A2
    data1 = get_auxillary_data(verts1)
    data2 = get_auxillary_data(verts2)
    Jac = newt_Jac(data1, data2)
    face_list = list_faces(data1, data2)
    ray_list = compute_discriminantal_normal_fan(data1, data2,Jac)
    Delta = cpt_discriminant_using_mixed_volumes(data1, data2, Jac, face_list, ray_list)
    Delta
end

# Given a planar polyhedral fan and a list of natural numbers, computes a polygon with prescribed normal fan.
# The lattice length of the edges are also prescribed.
#
# -------------  Input:
# ray_list              a List of Vectors
# ray_mults             a List of intergers
# -------------  Output:
# verts                a 2xn matrix

function pol_from_rays_and_mults(ray_list, ray_mults)
    verts = [0;0]
    for i in 1:lastindex(ray_list)-1
        ray = ray_list[i]
        mult = ray_mults[i]
        new_vert = verts[:,end] + mult*[ray[2],-ray[1]]
        verts = [verts new_vert]
    end
    m1 = minimum(verts[1,:])
    m2 = minimum(verts[2,:])
    verts[1,:] .-= m1
    verts[2,:] .-= m2
    verts
end

# Computes a polygon with prescribed normal fan, minimizing the total lattice length of the edges.
#
# -------------  Input:
# ray_list              a List of Vectors
# -------------  Output:
# pol                   a polyhedron

function get_minimal_polytope(ray_list)
    n_rays = length(ray_list)
    dualMat = hcat(ray_list...)
    dualMat = transpose(hcat(dualMat[2,:], -dualMat[1,:]))
    dualMat = MatrixSpace(ZZ, size(dualMat)...)(dualMat)

    
    -one(MatrixSpace(ZZ,n_rays,n_rays))
    A = vcat(dualMat,-dualMat,-one(MatrixSpace(ZZ,n_rays,n_rays)))
    b = vcat(zeros(Int64, 2*size(dualMat,1)),-ones(Int64,n_rays))
    cone = Polyhedron(A,b)
    LP = MixedIntegerLinearProgram(cone, ones(n_rays), convention = :min)
    _, ray_mults = solve_milp(LP)

    ray_mults = Int64.(lcm(denominator.(ray_mults)) * numerator.(Array(ray_mults)))

    verts = pol_from_rays_and_mults(ray_list, ray_mults)
    pol = convex_hull(transpose(verts))
    pol
end

# Compute a description of a polytope in terms of inequalities.
#
# -------------  Input:
# ray_list              a List of Vectors
# pol                   a polyhedron
# -------------  Output:
# weights                a list of integers

function get_inequalities(ray_list, pol)
    weights = Array{Int64}([])
    verts = [Array{Int64}(numerator.(vert)) for vert in vertices(pol)]
    for ray in ray_list
        weight = minimum([dot(ray, vert) for vert in verts])
        push!(weights, weight)
    end
    weights
end

# Compute the Newton polytope of g(f1,f2), where all polynomials are bivariate and generic withn prescribed Newton polytopes.
#
# -------------  Input:
# data1, data2             auxillary_data
# pol                      a polyhedron
# -------------  Output:
# P_result                 a polyhedron
function lifted_polytope(data1, data2, pol)
    verts = vertices(pol)
    vp1 = vertices(data1.pol)
    vp2 = vertices(data2.pol)
    MKsums = []
    for vert in verts
        MKsum = convex_hull( vert[1] .* vp1  ) + convex_hull( vert[2] .* vp2  )
        push!(MKsums, MKsum)
    end
    P_result = convex_hull(MKsums...)
    P_result
end

# Compute the Newton polytope of the Jacobian of a generic polynomial map f1, f2: C^2 -> C^2.
#
# -------------  Input:
# data1, data2             auxillary_data
# -------------  Output:
# Sigma                      a polyhedron
function newt_Jac(data1, data2)
    vert_list1 = data1.lpts
    vert_list2 = data2.lpts

    resulting_verts = []

    for vert1 in vert_list1
        for vert2 in vert_list2
            if det(hcat(vert1, vert2)) != 0
                push!(resulting_verts, vert1+vert2 - [1,1])
            end
        end
    end

    Sigma = Oscar.convex_hull( transpose(hcat(resulting_verts...)))
    Sigma
end



# Compute the Newton polytope of the discriminant of a generic polynomial map f1, f2: C^2 -> C^2.
# -------------  Input:
# data1, data2                      auxillary_data
# Sigma                             a polyhedron  
# face_list                         a List
# ray_list                          a List
# -------------  Output:
# Discriminant                      a polyhedron
function cpt_discriminant_using_mixed_volumes(data1, data2, Sigma, face_list, ray_list)
    correction = correction_term(face_list)
    weight_list = []
    test_pol = get_minimal_polytope(ray_list)
    A = transpose(hcat(ray_list...))
    b = get_inequalities(ray_list, test_pol)

    test_volume = taking_mixed_volume_with_discriminant(data1, data2, Sigma, test_pol, correction)
    for i in 1:length(ray_list)
        b_mod = copy(b)
        b_mod[i] += 1
        other_test_pol = Polyhedron(-A,-b_mod)
        #mixed_vol(correct_newt_D, test_pol) - mixed_vol(correct_newt_D, other_test_pol)
        #println("i = $i")
        weight = test_volume - taking_mixed_volume_with_discriminant(data1, data2, Sigma, other_test_pol, correction) 
        push!(weight_list, weight)
    end
    @assert(all(denominator.(weight_list) .== 1))
    weight_list = Int64.(numerator.(weight_list))

    verts = pol_from_rays_and_mults(ray_list, weight_list)
    Discriminant = convex_hull(transpose(verts))
    Discriminant
end

# Compute the mixed volume of a polyhedron with the Newton polytope of the discriminant of a generic polynomial map f1, f2: C^2 -> C^2.
# -------------  Input:
# data1, data2                      auxillary_data
# Sigma                             a polyhedron  
# pol                               a polyhedron
# correction                        a List
# -------------  Output:
# mvol                               a rational number
function taking_mixed_volume_with_discriminant(data1, data2, Sigma, pol, correction)
    verts = vertices(pol)
    transl = [ -minimum([vert[1] for vert in verts]) , -minimum([vert[2] for vert in verts]) ]
    pol_transl = pol+transl
    pol_lifted = lifted_polytope(data1, data2, pol_transl)
    mvol = volume(Sigma + pol_lifted) - volume(Sigma) - volume(pol_lifted)
    lpts = lattice_points(pol_transl)
    (a_max, b_max) = Int64.([ maximum([lpt[1] for lpt in lpts]) , maximum([lpt[2] for lpt in lpts]) ])
    mvol = (mvol - dot(correction, [a_max; b_max]))#[1]
    mvol
end



# Compute the normal fan of the Newton polytope of the discriminant of a generic polynomial map f1, f2: C^2 -> C^2.
# -------------  Input:
# data1, data2                      auxillary_data
# Sigma                             a polyhedron  
# -------------  Output:
# ray_list                          a List
function compute_discriminantal_normal_fan(data1, data2, Sigma)
    RR = TropicalSemiring(min)

    rays_of_Jacobian = rays(normal_fan(Sigma))

    Tf1_triv = data1.Support_Function
    Tf2_triv = data2.Support_Function

    new_rays = []
    for ray in rays_of_Jacobian
        vec = RR.(Array(ray))
        normal_ray = [Tf1_triv(vec...).data, Tf2_triv(vec...).data]
        if sum(abs.(normal_ray)) > 0
            normal_ray = normal_ray//gcd(normal_ray)
            normal_ray = (normal_ray.//lcm(denominator.(normal_ray)))
            push!(new_rays, Int.(normal_ray))
        end
    end
    push!(new_rays, [[1,0],[0,1]]...)

    ray_list = unique!(new_rays)
    sort!(ray_list, by = ray -> atan(ray[2], ray[1]))
end



function correction_term(face_list)
    correction = [0 0]
    edges = filter(face -> is_dicritical(face) && is_long(face) && !(is_origin(face)), face_list)
    for edge in edges
        γ1, γ2, ω = edge
        if !([0;0] in γ1)
            index = 1
            γ = γ1
        else
            index = 2
            γ = γ2
        end
        len = gcd(γ[2]-γ[1]) 
        correction[index] += len
    end
    correction
end

#  A data structure that stores useful information about a polyhedron
# -------------  Input:
# verts                   a 2xn matrix storing the vertices of the polyhedron.
# pol                     a polyhedron
# pol0                    the convex hull of pol with the origin
# lpts                    a list of lattice points
# m1                      the minimal x-coordinate of all vertices
# m2                      the minimal y-coordinate of all vertices
# ray_list                a list of rays in the normal fan
# no_origin               a boolean
# Support_Function        the support function of the polyhedron

struct auxillary_data
    verts
    pol
    pol0
    lpts
    m1
    m2
    ray_list
    no_origin
    Support_Function
end

# Computes useful information about a polyhedron
# -------------  Input:
# verts                 a 2xn matrix
# -------------  Output:
# data                  auxillary_data
function get_auxillary_data(verts)

    pol = convex_hull(transpose(verts))
    pol0 = convex_hull(transpose([verts [0;0]]))
    lpts = Vector{Vector{Int64}}(lattice_points(pol))

    S,(z1,z2) = PolynomialRing(ZZ,2)

    A = filter(vert -> !(vert[1]==0 && vert[2]==0), lpts)
    m1 = minimum([vert[1] for vert in A])
    m2 = minimum([vert[2] for vert in A])

    ray_list = []
    for ω in rays(normal_fan(pol0))
        ω_new = (1//gcd(ω)).*ω
        n_ray = Vector{Int64}([ω_new[1].num, ω_new[2].num])
        push!(ray_list, n_ray)
    end
    no_origin = all([!all(verts[:,j].== 0) for j in 1:size(verts,2)])


    RR = TropicalSemiring(min)
    V = transpose(verts)
    S,(x,y) = RR["x","y"]
    Support_Function = sum([   x^V[i, 1] * y^V[i, 2]  for i in 1:size(V, 1)  ])


    data = auxillary_data(verts, pol, pol0, lpts, m1, m2, ray_list, no_origin, Support_Function)
    data
end



# Lists all faces of a pair A = (A1, A2) of polygons. A face γ = (γ1, γ2) is represented by a triple consisting of:
# A list of vertices of γ1, a list of vertices of γ2 and a weight vector ω, revealing γ.
# -------------  Input:
# data1             auxillary_data
# data2             auxillary_data
# -------------  Output:
# face_list          list
function list_faces(data1, data2)
    verts1 = data1.verts
    verts1 = [verts1[:,x] for x in 1:size(verts1,2)]
    verts1 = [verts1; [[0;0]]]
    verts2 = data2.verts
    verts2 = [verts2[:,x] for x in 1:size(verts2,2)]
    verts2 = [verts2;[[0;0]]]


    ray_list1 = data1.ray_list
    ray_list2 = data2.ray_list
    ray_list = unique!([ray_list1; ray_list2])

    face_list = []
    for ω in ray_list
        γ1 = extremal_face(verts1, ω)
        γ2 = extremal_face(verts2, ω)
        if (length(γ1) > 1) | (length(γ2) > 1)
            push!(face_list,(γ1,γ2,ω))
        end
    end
    face_list
end



# Checks wether a face γ of A is dicritical
# -------------  Input:
# γ1                 a List
# -------------  Output:
# bool               boolean
function is_dicritical(face)
    (γ1, γ2, ω) = face
    bool = is_semi_origin(face) && ((ω[1] < 0) | (ω[2] < 0))
    bool
end

# Checks wether a face γ of A is long
# -------------  Input:
# γ1                 a List
# -------------  Output:
# bool               boolean
function is_long(face)
    (γ1, γ2, ω) = face
    length(γ1) > 1 && length(γ2) > 1
end

# Checks wether a face γ of A is origin
# -------------  Input:
# γ1                 a List
# -------------  Output:
# bool               boolean
function is_origin(face)
    (γ1, γ2, ω) = face
    Mγ1 = [[hcat(γ1...) [0; 0]]; ones(Int64,length(γ1)+1)']
    Mγ2 = [[hcat(γ2...) [0; 0]]; ones(Int64,length(γ2)+1)']

    is_semi_origin(face) && (rank([Mγ1 Mγ2]) ≤ 2)
end

# Checks wether a face γ of A is semi-origin
# -------------  Input:
# γ1                 a List
# -------------  Output:
# bool               boolean
function is_semi_origin(face)
    (γ1, γ2, ω) = face
    Mγ1 = [[hcat(γ1...) [0; 0]]; ones(Int64,length(γ1)+1)']
    Mγ2 = [[hcat(γ2...) [0; 0]]; ones(Int64,length(γ2)+1)']
    bool1 = (rank(Mγ1) == rank(Mγ1[:,1:end-1]))
    bool2 = (rank(Mγ2) == rank(Mγ2[:,1:end-1]))

    if length(γ1) == 1
        if γ1[1] != [0,0]
            false
        end
    end
    if length(γ2) == 1
        if γ2[1] != [0,0]
            false
        end
    end
    bool1 | bool2
end

# Checks wether a face γ of A is lower
# -------------  Input:
# γ1                 a List
# -------------  Output:
# bool               boolean
function is_lower(face)
    (γ1, γ2, ω) = face
    is_dicritical(face) && is_semi_origin(face) && (ω[2] > 0)
end

# Checks wether a face γ of A is upper
# -------------  Input:
# γ1                 a List
# -------------  Output:
# bool               boolean
function is_upper(face)
    (γ1, γ2, ω) = face
    is_dicritical(face) && is_semi_origin(face) && (ω[1] > 0)
end

# Give a weight vector ω, and a list of vertices, computes the vertices with minimal inner product.
# -------------  Input:
# verts                 a list
# -------------  Output:
# bovert_listl          a list
function extremal_face(verts, ω)
    min = minimum(transpose(ω)*hcat(verts...))
    indList = filter(i -> transpose(ω)*verts[i] == min, 1:length(verts))
    vert_list = verts[indList]
    vert_list
end

# Given two lists of vertices, computes the mixed volume of their respective convex hulls
# -------------  Input:
# verts1                a 2xn matrix
# verts2                a 2xn matrix
# -------------  Output:
# mvol                  an integer
function mixed_vol(verts1::Matrix, verts2::Matrix)
    mvol = Int64(Polymake.polytope.mixed_volume(convex_hull(transpose(verts1)).pm_polytope, convex_hull(transpose(verts2)).pm_polytope))
    mvol = Int64(numerator(mvol))
    mvol
end

# Given two polytopes, computes their mixed volume
# -------------  Input:
# pol1                a polyhedron
# pol2                a polyhedron
# -------------  Output:
# mvol                  an integer
function mixed_vol(pol1::Polyhedron, pol2::Polyhedron)
    mvol = (Polymake.polytope.mixed_volume(pol1.pm_polytope, pol2.pm_polytope))
    mvol = Int64(numerator(mvol))
    mvol
end

# For a generic polynomial map (f1, f2):C^2 -> C^2 computes the Newton polyhedra Gamma and Gamma'.
# -------------  Input:
# data1                 auxillary_data
# data2                 auxillary_data
# -------------  Output:
# Gamma                  a polyhedron
# Gamma'                 a polyhedron
function get_Gamma(data1,data2)
    face_list = list_faces(data1,data2)
    Gamma = get_Gamma(data1,data2,face_list,true)
    Gammap = get_Gamma(data1,data2,face_list,false)
    Gamma, Gammap
end

function get_Gamma(data1,data2,face_list,cpt_upper_faces::Bool)
    if cpt_upper_faces
        is_lower_upper = is_upper
    else
        is_lower_upper = is_lower
    end

    if !any(is_dicritical.(face_list).&&is_lower_upper.(face_list))
        return 0
    end
    lower_upper_long_edges = filter(face -> is_long(face) && is_lower_upper(face), face_list)
    if isempty(lower_upper_long_edges)
        return 0
    end
    edge = lower_upper_long_edges[1]
    γ1, γ2, ω = edge
    if !([0;0] in γ1)
        #return(convex_hull([γ1[2],γ1[1]]))
        #return
        l = abs(gcd(γ1[2]-γ1[1])) 
        return convex_hull([[0,0],[l,0]])
    elseif !([0;0] in γ2)
        l =  abs(gcd(γ2[2]-γ2[1])) 
        #return(convex_hull([γ2[2],γ2[1]]))
        return convex_hull([[0,0],[0,l]])
    end

    if data1.no_origin
        γ1 = γ1[1:end-1]
    end
    if data2.no_origin
        γ2 = γ2[1:end-1]
    end

    get_Γγ([γ1, γ2, ω])
end



# 
# -------------  Input:
# face = (γ1, γ2, w)                 a List
# γ1                                 a list of vectors    
# γ2                                 a list of vectors    
# w                                  a vector
# -------------  Output:
# Gammga_gamma                       a polyhedron
function get_Γγ(face)
    γ1, γ2, ω = face
    n1 = minimum( abs.(gcd.(γ1)) )
    m1 = maximum( abs.(gcd.(γ1)) )
    n2 = minimum( abs.(gcd.(γ2)) )
    m2 = maximum( abs.(gcd.(γ2)) )
    Gammga_gamma = convex_hull([[0;n1], [0;m1], [n2;0], [m2;0]])
    Gammga_gamma
end


# Computes the Newton number of a polygon
# -------------  Input:
# pol                 a polyhedron
# -------------  Output:
# N                  an integer
function get_Newton_nr(pol)
    x_min = -optimal_value(MixedIntegerLinearProgram(pol,[-1,0]))
    y_min = -optimal_value(MixedIntegerLinearProgram(pol,[0,-1]))
    vol_p = volume(pol)
    vol_p0 = volume(convex_hull(pol, convex_hull([x_min y_min])))
    v = vol_p0 - vol_p
    if v == 0
        return 0
    end
    verts = vertices(pol)
    a = minimum([vert[1] for vert in filter(vert->Int64(vert[2]) == y_min,verts)])
    b = minimum([vert[2] for vert in filter(vert->Int64(vert[1]) == x_min,verts)])
    N = 2 * v - a - b + 1
    N
end

# Computes the horizontal and vertical gap of a pair of polyhedra.
# -------------  Input:
# data1                 auxillary_data
# data2                 auxillary_data
# -------------  Output:
# mv                  an integer
# mh                  an integer
function get_m(data1,data2)
    m1A1 = data1.m1
    m2A1 = data1.m2
    m1A2 = data2.m1
    m2A2 = data2.m2
    (mv, mh) = (maximum([0, m1A1 + m1A2 - 1]), maximum([0, m2A1 + m2A2 - 1]))
    mv, mh
end



# Given two polygons, computes their polyhedral type.
# -------------  Input:
# verts1                 a list of vectors, the vertices a a polygon.
# verts2                 a list of vectors, the vertices a a polygon.
# -------------  Output:
# psi                    a list of integers

function get_polyhedral_type(verts1, verts2)
    data1 = get_auxillary_data(verts1)
    data2 = get_auxillary_data(verts2)
    psi = polyhedral_type(data1,data2)
    psi
end

# Given two polygons, computes their polyhedral type.
# -------------  Input:
# data1                 auxillary_data
# data2                 auxillary_data
# -------------  Output:
# psi                    a list of integers
function polyhedral_type(data1::auxillary_data,data2::auxillary_data)
    face_list = list_faces(data1, data2)
    Sigma = newt_Jac(data1, data2);
    ray_list = compute_discriminantal_normal_fan(data1, data2,Sigma)
    Delta = cpt_discriminant_using_mixed_volumes(data1, data2, Sigma, face_list, ray_list);
    Gamma, Gammap = get_Gamma(data1,data2);

    mv, mh = get_m(data1,data2)
    mc = 1
    if dim(Sigma) ≤ 0
        mc = 0
    end

    if Gamma == 0 || Gammap == 0
        Gamma_plus_Gammap = 0
        mvol_Gamma_Gammap = 0
    else
        Gamma_plus_Gammap = Gamma + Gammap
        mvol_Gamma_Gammap = mixed_vol(Gamma, Gammap)
    end

    if Gamma == 0
        NGamma = 0
        h0Gamma = 0
    else
        NGamma = get_Newton_nr(Gamma)
        h0Gamma =  number_of_interior_pts(Gamma)
    end

    if Gammap == 0
        NGammap = 0
        h0Gammap = 0
    else
        NGammap = get_Newton_nr(Gammap)
        h0Gammap =  number_of_interior_pts(Gammap)
    end
    if Gamma_plus_Gammap == 0
        NGamma_plus_Gammap = 0
        h0Gamma_plus_Gammap = 0
    else
        NGamma_plus_Gammap = get_Newton_nr(Gamma_plus_Gammap)
        h0Gamma_plus_Gammap =  number_of_interior_pts(Gamma_plus_Gammap)
    end


    (NDelta, NSigma) = get_Newton_nr.([Delta, Sigma])
    (h0Delta, h0Sigma) = number_of_interior_pts.([Delta, Sigma])


    delList = Int64.(ones(5))
    for i in 1:5
        if ([NDelta, NSigma, NGamma, NGammap, NGamma_plus_Gammap])[i] == 0
            delList[i] = 0
        end
    end
    (delDelta, delSigma, delGamma, delGammap, delGamma_plus_Gammap) = delList






    psi1 = mixed_vol(data1.pol0, data2.pol0)
    psi2 = h0Sigma
    psi3 = NSigma
    psi4 = mc+mh+mv
    psi5 = mc*mh+mc*mv+mh*mv
    psi6 = mc*mh*mv
    psi7 = h0Delta - h0Sigma + delDelta
    psi8 = mvol_Gamma_Gammap + NGamma_plus_Gammap - NGamma - NGammap
    psi9 = h0Gamma + h0Gammap + NGamma + NGammap
    psi10 = h0Gamma * h0Gammap + NGamma * NGammap
    psi11 = h0Gamma + h0Gammap + delGamma + delGammap
    psi12 = h0Gamma * h0Gammap + delGamma * delGammap

    psi = (psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8, psi9, psi10, psi11, psi12)
    psi = Int64.(numerator.(psi))
    psi
end

function number_of_interior_pts(pol)
    d =  dim(pol)
    if d ≤ 1
        return 0
    end
    nr = length(interior_lattice_points(pol))
    nr
end





function do_experiments(deg::Int)
    result_list = []
    pol_list = list_all_polygons(deg)
    conic_pols = filter(verts->is_conic(verts), pol_list)    
    data_list = get_auxillary_data.(conic_pols)

    for i in 1:length(data_list)
        u = ReentrantLock()
        Threads.@threads for j in i:length(data_list)
            (data1,data2) = data_list[[i,j]]
            new_psi = polyhedral_type(data1,data2)
            Threads.lock(u) do
                push!(result_list, ((data1,data2),new_psi))
            end
            #i = rand(1:length(data_list))
            #j = rand(1:length(data_list))
        end
    end
    result_list
end

deg = 2;
result_list = do_experiments(deg)

#lst = copy(result_list);
types = [res[2] for res in result_list];
#inds = sortperm(types);
#sorted_res = lst[inds];
appearing_types = unique(types);

#unique([type[(end-4):end] for type in appearing_types])


println("We consider the degree = $deg case.")
println("There are ", length(result_list), " many conical pairs")
println("These have ", length(appearing_types), " different topological types.")
