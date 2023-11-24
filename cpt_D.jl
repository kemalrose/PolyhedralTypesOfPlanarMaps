




using Oscar



#struct aux_data
#    V
#    verts
#    pol
#    pol0
#    lpts
#    f
#    ∇f
#    Tf
#    Tf_triv
#    mons
#    m1
#    m2
#    ray_list
#    no_origin
#end

#struct combined_data
#    Tf
#    TH
#    mons1
#    mons2
#    incidence1
#    incidence2
#    dict1
#    dict2
#    codim_1_cells
#    two_cells
#    verts
#    Jac
#    is_rel
#    is_lat
#    is_ess
#    is_adj
#    face_list
#end
#
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





#function combine_data(data1, data2, Vars)
#    RR = TropicalSemiring(min)
#    (x,y,z1,z2) = Vars
#    Tf1, Tf2 = data1.Tf, data2.Tf
#    trans_x, trans_y = RR.([rand(-100:100), rand(-100:100)])
#    Tf1_trans = Tf1(trans_x * x, trans_y * y)
#    Tf = Tf1_trans * Tf2
#
#    mons1, mons2 = monomials(Tf1_trans,Vars), data2.mons
#
#    f1 = rand(-10000:10000, length(data1.lpts))'*[z1^pt[1]*z2^pt[2] for pt in data1.lpts]
#    f2 = rand(-10000:10000, length(data2.lpts))'*[z1^pt[1]*z2^pt[2] for pt in data2.lpts]
#    Jac = det(Oscar.jacobi_matrix([f1,f2]))
#
#    TH, codim_1_cells, two_cells, verts = top_cells(Tf)
#    #if ! is_smooth(TH)
#    #    error("Have non transversal intersection of tropical hypersurfaces.")
#    #end
#
#    incidence1, dict1 = mon_incidence(mons1, two_cells, verts)
#    incidence2, dict2 = mon_incidence(mons2, two_cells, verts)
#
#    is_rel = is_relevant(mons1, mons2, incidence1, incidence2, dict1, dict2, two_cells, verts)
#    is_lat = is_lateral(mons1, mons2, incidence1, incidence2, dict1, dict2, two_cells, verts, is_rel)
#    is_ess, is_adj = is_essential(mons1, mons2, incidence1, incidence2, dict1, dict2, two_cells, verts, is_lat, is_rel, codim_1_cells)
#    face_list = list_faces(data1,data2)
#
#    #relev_inds = findall(is_rel)
#   combined_data(Tf, TH, mons1, mons2, incidence1, incidence2, dict1, dict2, codim_1_cells, two_cells, verts, Jac, is_rel, is_lat, is_ess, is_adj,face_list)
#end

#function trop_pols(V1, V2)
#    weight1 = weight_triang(V1)
#    weight2 = weight_triang(V2)
#    trans_x, trans_y = RR.([7, -11])
#    Tf1 = sum([   (x * trans_x)^V1[i, 1] * (y * trans_y)^V1[i, 2] * RR(weight1[i]) for i in 1:size(V1, 1)  ])
#    Tf2 = sum([   x^V2[i, 1] * y^V2[i, 2] * RR(weight2[i]) for i in 1:size(V2, 1)  ])
#    Tf = Tf1 * Tf2
#    TH = TropicalHypersurface(Tf)
##    cplx = TH.polyhedralComplex.pm_complex
#    #visualize(TH.polyhedralComplex)
#    Tf1, Tf2, Tf
#end





#function weight_triang(V)
#    dim_of_poly = rank([V ones(Int64,size(V,1))]) - 1
#    if dim_of_poly < 2
#        triang = [[collect(1:size(V,1))]]
#    else
#        triang = regular_triangulation(V)
#    end
#    if length(triang[1]) == 1
#        weight_vector = ones(fmpq, size(V,1))
#    else
#        incidence = IncidenceMatrix(triang[1])
#        SOP = SubdivisionOfPoints(V, incidence)
#        subdiv = SOP.pm_subdivision
#        weight_vector = subdiv.MIN_WEIGHTS
#    end
#    #secondary = secondary_cone(SOP::SubdivisionOfPoints)
#    #weights = sum(  [rays(secondary)[:]; lineality_space(secondary)[:]] ) 
#    #weights.//gcd(weights)
#    Vector{fmpq}(weight_vector).//gcd(weight_vector)
#end
#






#function trop_val(Tf, pt)
#    if pt[1] == 0
#        result = minimum(pt[2:3]'* Rational{Int64}.(Int64.(Tf.exps))[[2,1],:])
#    else
#        result = minimum(pt'* [Rational{Int64}.(Int64.(Tf.coeffs))';Rational{Int64}.(Tf.exps)[[2,1],:]])
#    end
#    result
#end




#function top_cells(Tf)
#    TH = TropicalHypersurface(Tf)
#    cplx = TH.polyhedralComplex.pm_complex
#    verts = convert(Matrix{Rational{Int64}}, cplx.VERTICES)
#
#
#    codim_1_cells = [findall(cplx.MAXIMAL_POLYTOPES[i, :]) for i in 1:size(cplx.MAXIMAL_POLYTOPES, 1)]
#
#    exps = [Int64.(Tf.exps[reverse(1:end),i]) for i in 1:size(Tf.exps, 2)]
#    weights = Int64.(Tf.coeffs)
#
#    two_cells = Vector{Int}[]
#    for i in 1:length(exps)
#        w, exp = weights[i], exps[i]
#        mon = x^exp[1] * y^exp[2] * RR(w)
#    
#        cell = []
#        for j in 1:size(verts, 1)
#            vert = verts[j,:]
#            if vert[1] == 1
#                val_Tf = Tf(vert[2], vert[3])
#                val_mon = mon(vert[2], vert[3])
#            else
#                val_Tf = sum([  prod(RR.(vert[2:3]).^ exps[l]) for l in 1:length(exps)])
#                val_mon = RR(vert[2])^exp[1] * RR(vert[3])^exp[2]
#            end
#            if val_Tf == val_mon
#                push!(cell, j)
#            end
#        end
#        if length(cell) > 1
#            push!(two_cells, cell)
#        end
#    end
#    TH, codim_1_cells, two_cells, verts
#end



#function mon_incidence(mons, cells, verts)
#    incidence = []
#    dict = Dict()
#
#    vals = zeros(Rational{Int64}, length(mons), size(verts,1))
#    for i in 1:length(mons)
#        mon = mons[i]
#        for j in 1:size(verts,1)
#            vert = verts[j,:]
#            vals[i,j] = trop_val(mon, vert)
#        end
#    end
#
#    TF_vals = [minimum(vals[:,j]) for j in 1:size(verts,1)]
#
#    for j in 1:length(cells)
#        cell = cells[j]
#        for i in 1:length(mons)
#            if(vals[i, cell] == TF_vals[cell])
#                push!(incidence, [i, j])
#                dict[j] = i
#                dict[cells[j]] = mons[i]
#            end
#        end
#    end
#    incidence, dict
#end
#
#function monomials(Tf, Vars)
#    (x,y,z1,z2) = Vars
#    [ x^Tf.exps[2, i] * y^Tf.exps[1, i] * Tf.coeffs[i]  for i in 1:length(Tf)  ] 
#end


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


#function new_ray(cell, cbd::combined_data)
#    (Tf,TH,mons1,mons2,incidence1,incidence2,dict1,dict2,two_cells,verts,Jac,is_rel,is_lat,is_ess,is_adj) = (cbd.Tf,cbd.TH,cbd.mons1,cbd.mons2,cbd.incidence1,cbd.incidence2,cbd.dict1,cbd.dict2,cbd.two_cells,cbd.verts,cbd.Jac,cbd.is_rel,cbd.is_lat,cbd.is_ess,cbd.is_adj)
#    two_cell_indices = filter(index -> issubset(cell, two_cells[index]), 1:length(two_cells))
#    j = two_cell_indices[1]
#    mon1 = mons1[dict1[j]]
#    mon2 = mons2[dict2[j]]
#
#    v1 = reverse(Int64.(mon1.exps))
#    v2 = reverse(Int64.(mon2.exps))
#    
#    inds = findall(i -> verts[cell[i],1] == 0, 1:length(cell))
#    res = verts[cell[inds], 2:end] * [v1 v2]
#    if rank(res) > 1
#        error("Applying Tf to cell gives something positive dimensional.")
#    end
#    #if rank(res) == 0
#    #    error("Applying Tf to cell gives back zero.")
#    #end
#    if rank(res) == 0
#        return zeros(Int64, 2)
#    end
#    index_nonzero_vector = findfirst(i -> sum(res[i,:].!=0)!=0, 1:size(res, 1))
#    ray = res[index_nonzero_vector,:]
#    ray = Vector{Int64}(ray/abs(gcd(ray)))
#end


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


#function get_rays(cbd::combined_data)
#
#    ray_list = Vector{Int64}.([[0; 1], [1; 0]])
#
#    for i in 1:length(cbd.codim_1_cells)
#        cell  = cbd.codim_1_cells[i]
#        if cbd.is_ess[i] 
#            ray = new_ray(cell, cbd)
#            push!(ray_list, ray)
#        end
#    end
#
#    for i in 1:length(cbd.codim_1_cells)
#        cell  = cbd.codim_1_cells[i]
#        if cbd.is_adj[i] 
#            ray = new_ray(cell, cbd)
#            push!(ray_list, ray)
#        end
#    end
#
#    for i in 1:length(cbd.two_cells)
#        cell = cbd.two_cells[i]
#        if cbd.is_rel[i] && !cbd.is_lat[i]
#            ray = new_ray(cell, cbd)
#            push!(ray_list, ray)
#        end
#    end
#    unique!(ray_list)
#    sort!(ray_list, by = ray -> atan(ray[2], ray[1]))
#end
#
#function get_ray_mults_old(cbd::combined_data, ray_list, data1::aux_data, data2::aux_data,Vars)    
#    (Tf,TH,mons1,mons2,incidence1,incidence2,dict1,dict2,two_cells,verts,Jac,is_rel,is_lat,is_ess,is_adj) = (cbd.Tf,cbd.TH,cbd.mons1,cbd.mons2,cbd.incidence1,cbd.incidence2,cbd.dict1,cbd.dict2,cbd.two_cells,cbd.verts,cbd.Jac,cbd.is_rel,cbd.is_lat,cbd.is_ess,cbd.is_adj);
#    (x,y,z1,z2) = Vars
#    
#    Lpts1 = data1.lpts
#    Lpts2 = data2.lpts
#    #(X, Y, U, W) = (@polyvar X[1:4])[1]
#
#    f1 = data1.f
#    f2 = data2.f
#
#    #f1 = randn(Float64, length(Lpts1))'*[X^pt[1]*Y^pt[2] for pt in Lpts1]
#    #f2 = randn(Float64, length(Lpts2))'*[X^pt[1]*Y^pt[2] for pt in Lpts2]
#    #Jac = differentiate(f1, X) * differentiate(f2, Y) - differentiate(f1, Y) * differentiate(f2, X)
#
#    correction = correction_term(data1, data2, cbd)
#
#
#    dual_ray_list = []
#    for ray in ray_list
#        push!(dual_ray_list, [-ray[2], ray[1]])
#        push!(dual_ray_list, [ray[2], -ray[1]])
#    end
#    sort!(dual_ray_list, by = ray -> atan(ray[2], ray[1]))
#    unique!(dual_ray_list)
#
#    vecs = [dual_ray_list[i]+dual_ray_list[i+1] for i in 1:length(dual_ray_list)÷2]
#    #vecs = dual_ray_list[1:end\div]
#    #vecs = [vecs ; [ 10*vec+[rand(-5:5),rand(-5:5)]  for vec in vecs]]
#
#    #vecs = [rand(-100:100, 2) for i in 1:round(length(ray_list))+1]
#    vecs = [vec.÷gcd(vec) for vec in vecs]
#
#    M_inters = zeros(Int64, length(ray_list), length(vecs))
#    is_pos = zeros(Bool, length(ray_list), length(vecs))
#    m_vols = zeros(Int64, length(vecs))
#
#    for i in 1:length(ray_list)
#        for j in 1:length(vecs)
#            vec = vecs[j]
#            ray = ray_list[i]
#            is_pos[i, j] = ray'*vec > 0
#            M_inters[i, j] = abs( ray[2] * vec[2] + ray[1] * vec[1] )
#            #abs(det([[-ray[2];ray[1]] vec]))
#        end
#    end
#
#    for i in 1:length(vecs)
#        vec = vecs[i]
#        b = [-min(0, vec[1]); -min(0, vec[2])]
#        a = vec+b
#
#        h = rand(-500:500)*z1^a[1]*z2^a[2]+rand(-500:500)*z1^b[1]*z2^b[2];
#        g = h(f1, f2)
#
#        #vol = mixed_volume([f1 - U, f2 - W, Jac, pol_C])
#        vol = Polymake.polytope.mixed_volume(newton_polytope(cbd.Jac).pm_polytope,newton_polytope(g).pm_polytope)
#        #vol = vol - 2*a[1]
#        vol = vol - sum(correction * [a b])
#        
#        #vol2 = Polymake.polytope.mixed_volume(newt_D.pm_polytope, newton_polytope(h).pm_polytope)
#        
#        #if vol != vol2
#        #    print(a[1],b[2])
#        #    @assert vol == vol2
#        #end
#        #R = parent(f1)
#        #I = ideal(R,[cbd.Jac,g]);
#        #I = (I:ideal(R,(z1*z2)^100));
#        #vdim(quo(R,I)[1])
#        m_vols[i] = vol 
#    end
#
#    M = transpose([M_inters.*is_pos M_inters.*(.!is_pos)])
#    M = MatrixSpace(ZZ, size(M)...)(M)
#    m = [m_vols; m_vols]
#    m = MatrixSpace(ZZ, length(m),1)(m)
#    if rank(M) < length(ray_list)
#        error("The incidence matrix is not of full rank!")
#    end
#    sol = Oscar.solve(M, m)
#    Matrix{Int64}(sol)[:]
#end
#
#
#function get_D(data1,data2,Jac, face_list)
#    #cbd = combine_data(data1,data2,Vars)
#    ray_list = get_rays_slim(data1, data2,Jac)
#    ray_mults =  get_ray_mults_fast(Jac, ray_list, data1, data2, face_list)  
#    verts = pol_from_rays_and_mults(ray_list, ray_mults)
#    convex_hull(transpose(verts))
#end
#
#function cpt_D_from_scratch(verts1, verts2)
#    #cbd = combine_data(data1,data2,Vars)
#    data1 = get_data_fast(verts1)[2]
#    data2 = get_data_fast(verts2)[2]
#    Jac = newt_Jac(data1, data2)
#    #cbd = combine_data(data1,data2,Vars);
#    face_list = list_faces(data1, data2)
#    ray_list = get_rays_slim(data1, data2,Jac)
#    ray_mults =  get_ray_mults_fast(Jac, ray_list, data1, data2, face_list)  
#    verts = pol_from_rays_and_mults(ray_list, ray_mults)
#    convex_hull(transpose(verts))
#end
#
#



function random_poly(newton_poly, R)
    allverts = Vector{Int64}.(lattice_points(newton_poly))
    coeffs = rand(-300000:300000, length(allverts))'
    coeffs*[R[1]^pt[1]*R[2]^pt[2] for pt in allverts]
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




#function get_ray_mults_fast(Jac, ray_list, data1, data2, face_list)
#
#    S, (z1,z2) = PolynomialRing(ZZ, 2)
#    correction = correction_term(data1, data2, face_list)
#
#    vecs = sample_from_cones(ray_list,5)
#
#    M = zeros(Int64, length(ray_list), 0)
#    m_vols = []
#
#    for j in 1:length(vecs)
#        vec = vecs[j]
#
#        M_lattice_inters, is_pos,  vol = cpt_next_linear_cond(Jac, ray_list, data1, data2, face_list, vec, correction, z1, z2)
#
#        M = hcat(M, [M_lattice_inters.*is_pos M_lattice_inters.*(.!is_pos)])
#
#        m_vols = [m_vols; ZZ.(numerator.([vol, vol]))]
#
#        if rank(M) == length(ray_list)
#            break
#        end
#    end
#
#
#
#    if rank(M) < length(ray_list)
#        verts1 = data1.verts
#        verts2 = data2.verts
#        println("Get error for:")
#        println("verts1 = $verts1")
#        println("verts2 = $verts2")
#        error("The incidence matrix is not of full rank!")
#    end
#
#
#    m = vcat(m_vols...)
#    M = transpose(MatrixSpace(ZZ, size(M)...)(M))
#    m = MatrixSpace(ZZ, length(m),1)(m)
#    sol = Oscar.solve(M, m)
#    Matrix{Int64}(sol)[:]
#end

#function cpt_next_linear_cond(Jac, ray_list, data1, data2, face_list, vec, correction, z1, z2)
#    (f1, f2) = (data1.f, data2.f)
#
#    is_pos = zeros(Bool, length(ray_list), 1)
#    M_lattice_inters = zeros(Int64, length(ray_list), 1)
#
#    for i in 1:length(ray_list)
#        ray = ray_list[i]
#        is_pos[i] = ray'*vec > 0
#        M_lattice_inters[i] = abs( ray[2] * vec[2] + ray[1] * vec[1] )
#    end
#
#
#    b = [-min(0, vec[1]); -min(0, vec[2])]
#    a = vec+b
#    h = z1^a[1]*z2^a[2]+z1^b[1]*z2^b[2];
#    g = h(f1, f2)
#    pol1 = Jac
#    pol2 = newton_polytope(g)
#    vol = volume(pol1 + pol2) - volume(pol1) - volume(pol2)
#
#    exps = hcat(collect(exponents(h))...)
#    a_max = maximum(exps[1,:])
#    b_max = maximum(exps[2,:])
#    vol = (vol - dot(correction, [a_max; b_max]))#[1]
#
#
#    #vol2 = Polymake.polytope.mixed_volume(correct_newt_D.pm_polytope, newton_polytope(h).pm_polytope)
#    #if vol != QQ(vol2)
#    #    print(a[1],b[2])
#    #    @assert vol == QQ(vol2)
#    #end
#
#
#    M_lattice_inters, is_pos, vol
#end
#
#function sample_from_cones(ray_list,k)
#    #dual_ray_list = []
#    #for ray in ray_list
#    #    push!(dual_ray_list, [-ray[2], ray[1]])
#    #    push!(dual_ray_list, [ray[2], -ray[1]])
#    #end
#    #sort!(dual_ray_list, by = ray -> atan(ray[2], ray[1]))
#    #unique!(dual_ray_list)
#    #for i in 1:length(dual_ray_list)÷2
##
#    #    for k in 1:5
#    #        vec = dual_ray_list[i]*(k-1) +dual_ray_list[i+1]*(5-k)
#    #    push!(vecs,vec)
##
#    #    end
#    #end
#
#    vecs = [ Array{Int}(vec) for vec in lattice_points(k*cube(2))]
#    filter!(vec-> sum(abs.(vec))>0,vecs)
#    for i in 1:length(vecs)
#        vec = vecs[i]
#        vec = vec .÷ gcd(vec)
#        if vec[1] < 0
#            vec = -1 .*vec
#        end
#        vecs[i] = vec
#    end
#    vecs = shuffle!(unique!(vecs))
#    vecs
#end



#vecs = filter(vec -> all([dot(vec, ray) != 0 for ray in ray_list]),vecs)
#shuffle!(vecs)



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


#function mixed_volume_via_pullback_old(data1, data2, Jac, pol, correction)
#
#    verts = vertices(pol)
#    transl = [ -minimum([vert[1] for vert in verts]) , -minimum([vert[2] for vert in verts]) ]
#    pol_transl = pol+transl
#
#
#    S,(z1,z2) = PolynomialRing(ZZ,2)
#
#    lpts = lattice_points(pol)
#    transl = [ -minimum([lpts[1] for lpts in lpts]) , -minimum([lpts[2] for lpts in lpts]) ]
#    lpts = [lpt + transl for lpt in lpts]
#    h = rand(1:1, length(lpts))'*[z1^pt[1]*z2^pt[2] for pt in lpts]
#    (f1, f2) = (data1.f, data2.f)
#    g = h(f1, f2)
#    pol2 = newton_polytope(g)
#    
#
#    pol1 = Jac
#    pol2_new = lifted_polytope(data1, data2, pol_transl)
#
#
#    vol = volume(pol1 + pol2) - volume(pol1) - volume(pol2)
#
#    vol_new = volume(pol1 + pol2_new) - volume(pol1) - volume(pol2_new)
#    vol = vol_new
#    if !(vol == vol_new)
#        verts_pol = vertices(pol)
#        println("verts = $verts_pol")
#    end
#    @assert vol == vol_new
#
#
#    exps = hcat(collect(exponents(h))...)
#    a_max = maximum(exps[1,:])
#    b_max = maximum(exps[2,:])
#    vol = (vol - dot(correction, [a_max; b_max]))#[1]
#
#    vol2 = Polymake.polytope.mixed_volume(correct_newt_D.pm_polytope, newton_polytope(h).pm_polytope)
#    
#    if !(vol == QQ(vol2))
#        verts_pol = vertices(pol_transl)
#        println("verts = $verts_pol")
#    end
#    @assert(vol == QQ(vol2))
#
#    vol
#end


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



function get_inequalities(ray_list, pol)
    weights = Array{Int64}([])
    verts = [Array{Int64}(numerator.(vert)) for vert in vertices(pol)]
    for ray in ray_list
        weight = minimum([dot(ray, vert) for vert in verts])
        push!(weights, weight)
    end
    weights
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

function get_D(A1, A2)
    verts1, verts2 = A1, A2
    data1 = get_data_fast(verts1)[2]
    data2 = get_data_fast(verts2)[2]
    Jac = newt_Jac(data1, data2)
    face_list = list_faces(data1, data2)
    ray_list = get_rays_slim(data1, data2,Jac)
    D = cpt_discriminant_chow(data1, data2, Jac, face_list, ray_list)
    D
end


#function get_discriminant(verts1, verts2)
#    bool, data1 = get_data_fast(verts1)
#    bool, data2 = get_data_fast(verts2)
#    Jac = newt_Jac(data1, data2)
#    face_list = list_faces(data1, data2)
#    ray_list = get_rays_slim(data1, data2,Jac)
#    cpt_discriminant_chow(data1, data2, Jac, face_list, ray_list)
#end



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