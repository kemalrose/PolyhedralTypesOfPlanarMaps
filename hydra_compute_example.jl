

#using Pkg
#Pkg.add("Oscar")
using Oscar, LinearAlgebra, Random


function list_all_polytopes(deg)
    S = [[i,j] for i in 0:deg for j in 0:deg]
    filter!(s-> sum(s) ≤ deg && sum(s) > 0, S)
    S = sort(S)
    polytopeList = []
    newPolytopes = []
    for s1 in S
        for s2 in S
            if s1 < s2
                push!(newPolytopes, hcat(s1, s2))
            end
        end
    end
    while length(newPolytopes) > 0
        polytopeList = [polytopeList; newPolytopes]
        newPolytopes = extendPolytopes(newPolytopes, S)
    end
    polytopeList
end



#for i in 1:nverts
#    s = pol[:,(i-1)%nverts+1]
#    p = pol[:,((i-1)+1)%nverts+1]
#    q = pol[:,((i-1)+2)%nverts+1]
#    push!(res,isAdmissible(s,p,q))
#end

#function isAdmissible(P, v)
#    s = P[:,1]
#    p = P[:,end-1]
#    q = P[:,end]
#    if size(P,2) == 2
#        return isAdmissible(s, q, v)
#    else
#
#        return isAdmissible(p, q, v) && isAdmissible(s, p, v) && s < v
#    end
#end


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



function isAdmissible(p, q, v)
    qn = q - p
    vn = v - p
    qn[2] * vn[1] - qn[1] * vn[2] > 0
end


function extendPolytopes(pols, S)
    result = []
    for pol in pols
        for s in S
            if isAdmissible(pol, s)
                push!(result, hcat(pol,s))
            end
        end
    end
    result
end


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

function is_conic(lpts)
    M = zeros(Int64,5,size(lpts,1))
    for j in 1:size(lpts,1)
        s1, s2 = lpts[j]
        M[1,j] = s1
        M[2,j] = s2
        M[3,j] = s1^2
        M[4,j] = s1*s2
        M[5,j] = s2^2
    end
    rank(M) >= 5
end

#function is_conic_old_pol(verts, lpts)
#    if size(verts,2) < 3
#        return false
#    end
#    if size(verts,2) > 4
#        return true
#    end
#    return is_conic_lpts(lpts)
#end
#
#function is_conic_lpts(lpts)
#    n_pts = lastindex(lpts)
#    index_list = [1]
#
#    while length(index_list) > 0
#        if length(index_list) > 4
#            return true
#        end
#        bool_found_new_index, extending_index = can_extend(index_list, lpts)
#        if bool_found_new_index
#            push!(index_list, extending_index)
#        else
#            index_list[end]+=1
#            if index_list[end] >= n_pts
#                pop!(index_list)
#                if length(index_list) > 0
#                    index_list[end]+=1
#                end
#            end
#        end
#    end
#    false
#end

function is_ok_to_add(index_list, lpts)
    for i in 1:lastindex(index_list)-1
        for j in i+1:lastindex(index_list)-1
            if areColinear(lpts[[i, j, index_list[end]]]...)
                return false
            end
        end
    end
    true
end

function can_extend(index_list, lpts)
    n_pts = lastindex(lpts)
    bool_found_new_index = false
    extending_index = nothing
    new_indices = collect(index_list[end]+1:n_pts)
    for index in new_indices
        if is_ok_to_add([index_list;index], lpts)
            extending_index = index
            bool_found_new_index = true
            break
        end
    end
    bool_found_new_index, extending_index
end




function areColinear(v1, v2, v3)
    rank(hcat(v1-v3, v2-v3)) < 2
end


function degree(pol)
    maximum([sum(pol[:,i]) for i in axes(pol,2)])
end



function get_delta(A1, A2)
    verts1, verts2 = A1, A2
    data1 = get_data_fast(verts1)[2]
    data2 = get_data_fast(verts2)[2]
    Jac = newt_Jac(data1, data2)
    face_list = list_faces(data1, data2)
    ray_list = get_rays_slim(data1, data2,Jac)
    D = cpt_discriminant_chow(data1, data2, Jac, face_list, ray_list)
    D
end


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


function get_newton_polytope(ray_list)
    n_rays = length(ray_list)
    dualMat = hcat(ray_list...)
    dualMat = transpose(hcat(dualMat[2,:], -dualMat[1,:]))
    dualMat = MatrixSpace(ZZ, size(dualMat)...)(dualMat)
    #Mat = transpose(nullspace(dualMat)[2])
    #ker = positive_hull(hcat(Mat, -Mat))
    -one(MatrixSpace(ZZ,n_rays,n_rays))
    A = vcat(dualMat,-dualMat,-one(MatrixSpace(ZZ,n_rays,n_rays)))
    b = vcat(zeros(Int64, 2*size(dualMat,1)),-ones(Int64,n_rays))
    cone = Polyhedron(A,b)
    LP = MixedIntegerLinearProgram(cone, ones(n_rays), convention = :min)
    m, ray_mults = solve_milp(LP)

    ray_mults = Int64.(lcm(denominator.(ray_mults)) * numerator.(Array(ray_mults)))

    #ray_mults = Array{Int64}(ray_mults)
    #kerplus = Polyhedron(-I[1:8, 1:8], -ones(Int64,8))
    #cone = intersect(ker, positive_hull(I[1:8, 1:8]))
    #intersect(kerplus, cone)

    #ray_mults
    verts = pol_from_rays_and_mults(ray_list, ray_mults)
    convex_hull(transpose(verts))
end



function get_inequalities(ray_list, pol)
    weights = Array{Int64}([])
    verts = [Array{Int64}(numerator.(vert)) for vert in vertices(pol)]
    for ray in ray_list
        weight = minimum([dot(ray, vert) for vert in verts])
        push!(weights, weight)
    end
    weights
end


function lifted_polytope(data1, data2, pol)
    verts = vertices(pol)
    vp1 = vertices(data1.pol)
    vp2 = vertices(data2.pol)
    MKsums = []
    for vert in verts
        MKsum = convex_hull( vert[1] .* vp1  ) + convex_hull( vert[2] .* vp2  )
        push!(MKsums, MKsum)
    end
    result = convex_hull(MKsums...)
    result
end


function mixed_volume_via_pullback(data1, data2, Jac, pol, correction)
    
    verts = vertices(pol)
    transl = [ -minimum([vert[1] for vert in verts]) , -minimum([vert[2] for vert in verts]) ]
    pol_transl = pol+transl
    pol_lifted = lifted_polytope(data1, data2, pol_transl)
    vol = volume(Jac + pol_lifted) - volume(Jac) - volume(pol_lifted)
    lpts = lattice_points(pol_transl)
    (a_max, b_max) = Int64.([ maximum([lpt[1] for lpt in lpts]) , maximum([lpt[2] for lpt in lpts]) ])
    vol = (vol - dot(correction, [a_max; b_max]))#[1]
    vol
end

function cpt_discriminant_chow(data1, data2, Jac, face_list, ray_list)
    correction = correction_term(data1, data2, face_list)
    weight_list = []
    test_pol = get_newton_polytope(ray_list)
    A = transpose(hcat(ray_list...))
    b = get_inequalities(ray_list, test_pol)

    test_volume = mixed_volume_via_pullback(data1, data2, Jac, test_pol, correction)
    for i in 1:length(ray_list)
        b_mod = copy(b)
        b_mod[i] += 1
        other_test_pol = Polyhedron(-A,-b_mod)
        #mixed_vol(correct_newt_D, test_pol) - mixed_vol(correct_newt_D, other_test_pol)
        #println("i = $i")
        weight = test_volume - mixed_volume_via_pullback(data1, data2, Jac, other_test_pol, correction) 
        push!(weight_list, weight)
    end
    @assert(all(denominator.(weight_list) .== 1))
    weight_list = Int64.(numerator.(weight_list))

    verts = pol_from_rays_and_mults(ray_list, weight_list)
    convex_hull(transpose(verts))
end


function random_poly(newton_poly, R)
    allverts = Vector{Int64}.(lattice_points(newton_poly))
    coeffs = rand(-300000:300000, length(allverts))'
    coeffs*[R[1]^pt[1]*R[2]^pt[2] for pt in allverts]
end




function get_D_elim(verts1, verts2)
    p = Primes.prime(100000)
    FF,_ = FiniteField(p,1,"")

    R,(z1,z2,y1,y2) = PolynomialRing(FF,["z1","z2","y1","y2"]);
    f1 = random_poly( convex_hull(verts1') ,R);
    f2 = random_poly(convex_hull(verts2'),R);
    Jac = det(jacobi_matrix([f1,f2])[1:2,:]);
    I = ideal(R,[Jac, f1-y1, f2-y2]);
    D_elim = eliminate(I,[z1,z2])[1];
    convex_hull([vert[3:4] for vert in vertices(newton_polytope(D_elim))])
end


function newt_Jac(data1, data2)
#    vert_list1 = [lpts1[:,i] for i in axes(lpts1)[2]]
#    vert_list2 = [lpts2[:,i] for i in axes(lpts2)[2]]
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

    Oscar.convex_hull( transpose(hcat(resulting_verts...)))
end



function get_rays_slim(data1, data2, Jac)
    RR = TropicalSemiring(min)
    NJac = Jac
    rays_C = rays(normal_fan(NJac))

    Tf1_triv = data1.Tf_triv
    Tf2_triv = data2.Tf_triv

    new_rays = []
    for ray in rays_C
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



function correction_term(data1, data2, face_list)
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



function is_relevant(mons1, mons2, incidence1, incidence2, dict1, dict2, two_cells, verts)
    bools = []

    for j in 1:length(two_cells)
        cell = two_cells[j]

        mon1 = mons1[dict1[j]]
        mon2 = mons2[dict2[j]]

        v1 = reverse(Int64.(mon1.exps))
        v2 = reverse(Int64.(mon2.exps))

        bool = false

        inds = findall(i -> verts[cell[i],1] == 0, 1:length(cell))
        res = verts[cell[inds], 2:end] * [v1 v2]
        if rank(hcat(v1, v2)) < 2 && !all(res.==0)               
            bool = true
        end
        push!(bools, bool)
    end

    bools
end

function is_bounded(cell, verts)
    sum(verts[cell,1].==0) == 0
end

function is_lateral(mons1, mons2, incidence1, incidence2, dict1, dict2, two_cells, verts, is_rel)
    bools = []

    all_exp_vectors = [  Int64.(mon.exps)  for mon in [mons1; mons2]]


    for j in 1:length(two_cells)
        cell = two_cells[j]

        mon1 = mons1[dict1[j]]
        mon2 = mons2[dict2[j]]

        v1 = Int64.(mon1.exps)
        v2 = Int64.(mon2.exps)

        bool = false
        if is_rel[j]
            v = v1 + v2
            inner_products = [v[2]*exp_vec[1] - v[1]*exp_vec[2] for exp_vec in all_exp_vectors]
            bool = all( inner_products.≥0) || all( inner_products.≤0)
        end

        push!(bools, bool)
    end
    bools
end

function is_essential(mons1, mons2, incidence1, incidence2, dict1, dict2, two_cells, verts, is_lat, is_rel, codim_1_cells)
    n_1_cells = length(codim_1_cells)
    n_2_cells = length(two_cells)
    bools_ess = []
    bools_adj = []

    for i in 1:n_1_cells
        bool_ess = false
        bool_adj = false
        line = codim_1_cells[i]
        if !is_bounded(line, verts)
            two_cell_indices = filter(index -> issubset(line, two_cells[index]), 1:n_2_cells)
            if (sum(is_rel[two_cell_indices]) == 0) 
                bool_adj = true
            end

            if (sum(is_rel[two_cell_indices]) == 1) && (sum(is_lat[two_cell_indices]) > 0)

                two_mons = [dict1[two_cell_indices[1]], dict1[two_cell_indices[2]]]
                two_mons = mons1[two_mons]
                if length(unique(two_mons)) != 2
                    two_mons = [dict2[two_cell_indices[1]], dict2[two_cell_indices[2]]]
                    two_mons = mons2[two_mons]
                end

                if is_lat[two_cell_indices[2]]
                    reverse!(two_mons)
                end
                v1, v2 = [Int64.(mon.exps) for mon in two_mons]
                w1 = (v2-v1)#//gcd(v2-v1)
                w2 = v1//gcd(v1)
                if abs(det(hcat(w1, w2))) > 1
                    bool_ess = true
                end
            end
        end
        push!(bools_ess, bool_ess)
        push!(bools_adj, bool_adj)


#        for j in 1:n_2_cells
#            bool = false
#            if is_lat[j]
#
#            end
#            boolean_array[i,j] = bool
#        end
    end
    bools_ess, bools_adj
end


struct data_fast
    f
    verts
    pol
    pol0
    lpts
    m1
    m2
    ray_list
    no_origin
    Tf_triv
end


function get_data_fast(verts)

    pol = convex_hull(transpose(verts))
    pol0 = convex_hull(transpose([verts [0;0]]))
    lpts = Vector{Vector{Int64}}(lattice_points(pol))
    bool = is_conic(lpts)

    S,(z1,z2) = PolynomialRing(ZZ,2)
    f = rand(1:1, length(lpts))'*[z1^pt[1]*z2^pt[2] for pt in lpts]

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
    Tf_triv = sum([   x^V[i, 1] * y^V[i, 2]  for i in 1:size(V, 1)  ])


    data = data_fast(f,verts, pol, pol0, lpts, m1, m2, ray_list, no_origin, Tf_triv)
    bool, data
end


function is_dicritical(face)
    (γ1, γ2, ω) = face
    is_semi_origin(face) && ((ω[1] < 0) | (ω[2] < 0))
end

function is_long(face)
    (γ1, γ2, ω) = face
    length(γ1) > 1 && length(γ2) > 1
end


function is_origin(face)
    (γ1, γ2, ω) = face
    Mγ1 = [[hcat(γ1...) [0; 0]]; ones(Int64,length(γ1)+1)']
    Mγ2 = [[hcat(γ2...) [0; 0]]; ones(Int64,length(γ2)+1)']

    is_semi_origin(face) && (rank([Mγ1 Mγ2]) ≤ 2)
end


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

function is_lower(face)
    (γ1, γ2, ω) = face
    is_dicritical(face) && is_semi_origin(face) && (ω[2] < 0)
end

function is_upper(face)
    (γ1, γ2, ω) = face
    is_dicritical(face) && is_semi_origin(face) && (ω[1] < 0)
end


function extremal_face(verts, ω)
    min = minimum(transpose(ω)*hcat(verts...))
    indList = filter(i -> transpose(ω)*verts[i] == min, 1:length(verts))
    verts[indList]
end

function restrictfct(allverts, coeffsf, ω)
    min = minimum(transpose(ω)*hcat(allverts...) )
    indList = filter(i -> transpose(ω)*allverts[i] == min, 1:length(allverts))
    fω = transpose(coeffsf[indList]) * [z1^pt[1]*z2^pt[2] for pt in allverts[indList]]
    fω, min
end

function mixed_vol(verts1, verts2)
    Polymake.polytope.mixed_volume(convex_hull(verts1).pm_polytope, convex_hull(verts2).pm_polytope)
end

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

function get_LU(data1,data2)
    face_list = list_faces(data1,data2)
    L = get_LU(data1,data2,face_list,false)
    U = get_LU(data1,data2,face_list,true)
    L, U
end

function get_LU(data1,data2,face_list,cpt_upper_faces::Bool)
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
        return -1
    end
    edge = lower_upper_long_edges[1]
    γ1, γ2, ω = edge
    if !([0;0] in γ1)
        return - gcd(γ1[2]-γ1[1]) - 1
    elseif !([0;0] in γ2)
        return - gcd(γ2[2]-γ2[1]) - 1
    end

    if data1.no_origin
        γ1 = γ1[1:end-1]
    end
    if data2.no_origin
        γ2 = γ2[1:end-1]
    end

    Γγ = get_Γγ([γ1, γ2, ω])
    σ1 = 0
    N = get_Newton_nr(Γγ)
    if N == 1
        σ1 += 1
    end
    length(interior_lattice_points(Γγ)) + σ1 + 1
end




function get_Γγ(face)
    γ1, γ2, ω = face
    n1 = minimum( abs.(gcd.(γ1)) )
    m1 = maximum( abs.(gcd.(γ1)) )
    n2 = minimum( abs.(gcd.(γ2)) )
    m2 = maximum( abs.(gcd.(γ2)) )
    convex_hull([[n1;0], [m1;0], [0;n2], [0;m2]])
end



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
    2 * v - a - b + 1
end

#function random_poly(newton_poly, R)
#    allverts = Vector{Int64}.(lattice_points(newton_poly))
#    coeffs = rand(-50:50, length(allverts))'
#    coeffs*[R[1]^pt[1]*R[2]^pt[2] for pt in allverts]
#end


function get_m(data1,data2)
    m1A1 = data1.m1
    m2A1 = data1.m2
    m1A2 = data2.m1
    m2A2 = data2.m2
    maximum([0, m1A1 + m1A2 - 1]), maximum([0, m2A1 + m2A2 - 1])
end



function polyhedral_type(data1::data_fast,data2::data_fast)
    #ray_list = get_rays_slim(data1, data2)
    #ray_list =  get_rays(cbd)

    Jac = newt_Jac(data1, data2)
    face_list = list_faces(data1, data2)
    ray_list = get_rays_slim(data1, data2,Jac)
    D = cpt_discriminant_chow(data1, data2, Jac, face_list, ray_list)

    mv, mh = get_m(data1,data2)
    N0 = get_Newton_nr(D);
    σ1 = 0
    σ2 = 0
    if N0 == 1
        σ1 += 1    
    elseif N0 == 2
        σ2 += 1
    end
    LA, UA =  get_LU(data1,data2);
    ψ1 = Polymake.polytope.mixed_volume(data1.pol0.pm_polytope, data2.pol0.pm_polytope) 
    ψ2 = length(interior_lattice_points(Jac))
    ψ3 = mv+mh
    ψ4 = mv*mh
    ψ5 = length(interior_lattice_points(D)) - length(interior_lattice_points(Jac)) + σ1 + σ2;
    ψ6 = LA + UA
    ψ7 = LA * UA
    ψ = (ψ1, ψ2, ψ3, ψ4, ψ5, ψ6, ψ7)
    Int64.(numerator.(ψ))
end



#function do_experiments(deg::Int)
#    psi_list = []
#    all_pols = listAllPolytopes(deg)
#    data_list = []
#    for verts in all_pols
#        bool, data = get_data_fast(verts)
#        if bool
#            push!(data_list, data)
#        end
#    end
#
#    #for index in 1:200
#
#    
#    for i in 1:length(data_list)
#        Threads.@threads for j in i:length(data_list)
#            #i = rand(1:length(data_list))
#            #j = rand(1:length(data_list))
#            (data1,data2) = data_list[[i,j]]
#            new_psi = polyhedral_type(data1,data2)
#            push!(psi_list, new_psi)
#        end
#    end
#    psi_list
#end

function get_polyhedral_type(A1,A2)
    verts1, verts2 = A1, A2
    data1 = get_data_fast(verts1)[2]
    data2 = get_data_fast(verts2)[2]
    polyhedral_type(data1,data2)
end


function list_all_conical_pols(deg::Int)
    all_pols = list_all_polytopes(deg)
    data_list = []
    for verts in all_pols
        bool, data = get_data_fast(verts)
        if bool
            push!(data_list, data)
        end
    end
    data_list
end


function do_experiments(deg::Int)
    result_list = []
    data_list = list_all_conical_pols(deg)

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



result_list = do_experiments(2);

lst = copy(result_list);
types = [res[2] for res in result_list];
inds = sortperm(types);
sorted_res = lst[inds];
appearing_types = unique(types);

println("There are ", length(types), " many conical pairs")
println("These have ", length(appearing_types), " different topological types.")

#for res in sorted_res
#    ((data1, data2),psi) = res
#    A1 = data1.verts
#    A2 = data2.verts
#    println("")
#    println("A1 = $A1")
#    println("A2 = $A2")
#    println("psi = $psi")
#end


#result_list = do_experiments(4)
#data_list = list_all_conical_pols(3)
#n_types = 0
#for i in 1:length(data_list)
#    for j in i:length(data_list)
#        n_types+=1
#    end
#end
#n_types