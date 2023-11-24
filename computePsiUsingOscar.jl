
using Oscar

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



function do_experiments(deg::Int)
    psi_list = []
    all_pols = listAllPolytopes(deg)
    data_list = []
    for verts in all_pols
        bool, data = get_data_fast(verts)
        if bool
            push!(data_list, data)
        end
    end

    #for index in 1:200

    
    for i in 1:length(data_list)
        Threads.@threads for j in i:length(data_list)
            #i = rand(1:length(data_list))
            #j = rand(1:length(data_list))
            (data1,data2) = data_list[[i,j]]
            new_psi = polyhedral_type(data1,data2)
            push!(psi_list, new_psi)
        end
    end
    psi_list
end
function get_polyhedral_type(A1,A2)
    verts1, verts2 = A1, A2
    data1 = get_data_fast(verts1)[2]
    data2 = get_data_fast(verts2)[2]
    polyhedral_type(data1,data2)
end
