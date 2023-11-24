


using Oscar, HomotopyContinuation

function trop_pols(V1, V2)
    weight1 = weight_triang(V1)
    weight2 = weight_triang(V2)
    trans_x, trans_y = RR.([7, -11])
    Tf1 = sum([   (x * trans_x)^V1[i, 1] * (y * trans_y)^V1[i, 2] * RR(weight1[i]) for i in 1:size(V1, 1)  ])
    Tf2 = sum([   x^V2[i, 1] * y^V2[i, 2] * RR(weight2[i]) for i in 1:size(V2, 1)  ])
    Tf = Tf1 * Tf2
    TH = TropicalHypersurface(Tf)
#    cplx = TH.polyhedralComplex.pm_complex
    #visualize(TH.polyhedralComplex)
    Tf1, Tf2, Tf
end


function weight_triang(V)
    dim_of_poly = rank([V ones(Int64,size(V,1))]) - 1
    if dim_of_poly < 2
        triang = [[collect(1:size(V,1))]]
    else
        triang = regular_triangulation(V)
    end
    if length(triang[1]) == 1
        weight_vector = ones(fmpq, size(V,1))
    else
        incidence = IncidenceMatrix(triang[1])
        SOP = SubdivisionOfPoints(V, incidence)
        subdiv = SOP.pm_subdivision
        weight_vector = subdiv.MIN_WEIGHTS
    end
    #secondary = secondary_cone(SOP::SubdivisionOfPoints)
    #weights = sum(  [rays(secondary)[:]; lineality_space(secondary)[:]] ) 
    #weights.//gcd(weights)
    Vector{fmpq}(weight_vector).//gcd(weight_vector)
end







#TH1 = TropicalHypersurface(Tf1)
#TH2 = TropicalHypersurface(Tf2)
#subdiv = dual_subdivision(TH)
#
#cplx1 = TH1.polyhedralComplex.pm_complex
#cplx2 = TH2.polyhedralComplex.pm_complex
#one_cells = [    findall(cplx.MAXIMAL_POLYTOPES[i, :]) for i in 1:size(cplx.MAXIMAL_POLYTOPES, 1) ]



function trop_val(Tf, pt)
    if pt[1] == 0
        result = minimum(pt[2:3]'* Rational{Int64}.(Int64.(Tf.exps))[[2,1],:])
    else
        result = minimum(pt'* [Rational{Int64}.(Int64.(Tf.coeffs))';Rational{Int64}.(Tf.exps)[[2,1],:]])
    end
    result
end




function top_cells(Tf)
    TH = TropicalHypersurface(Tf)
    cplx = TH.polyhedralComplex.pm_complex
    verts = convert(Matrix{Rational{Int64}}, cplx.VERTICES)


    codim_1_cells = [findall(cplx.MAXIMAL_POLYTOPES[i, :]) for i in 1:size(cplx.MAXIMAL_POLYTOPES, 1)]

    exps = [Int64.(Tf.exps[reverse(1:end),i]) for i in 1:size(Tf.exps, 2)]
    weights = Int64.(Tf.coeffs)

    two_cells = Vector{Int}[]
    for i in 1:length(exps)
        w, exp = weights[i], exps[i]
        mon = x^exp[1] * y^exp[2] * RR(w)
    
        cell = []
        for j in 1:size(verts, 1)
            vert = verts[j,:]
            if vert[1] == 1
                val_Tf = Tf(vert[2], vert[3])
                val_mon = mon(vert[2], vert[3])
            else
                val_Tf = sum([  prod(RR.(vert[2:3]).^ exps[l]) for l in 1:length(exps)])
                val_mon = RR(vert[2])^exp[1] * RR(vert[3])^exp[2]
            end
            if val_Tf == val_mon
                push!(cell, j)
            end
        end
        if length(cell) > 1
            push!(two_cells, cell)
        end
    end
    TH, codim_1_cells, two_cells, verts
end



function mon_incidence(mons, cells, verts)
    incidence = []
    dict = Dict()

    vals = zeros(Rational{Int64}, length(mons), size(verts,1))
    for i in 1:length(mons)
        mon = mons[i]
        for j in 1:size(verts,1)
            vert = verts[j,:]
            vals[i,j] = trop_val(mon, vert)
        end
    end

    TF_vals = [minimum(vals[:,j]) for j in 1:size(verts,1)]

    for j in 1:length(cells)
        cell = cells[j]
        for i in 1:length(mons)
            if(vals[i, cell] == TF_vals[cell])
                push!(incidence, [i, j])
                dict[j] = i
                dict[cells[j]] = mons[i]
            end
        end
    end
    incidence, dict
end

function monomials(Tf)
    [ x^Tf.exps[2, i] * y^Tf.exps[1, i] * Tf.coeffs[i]  for i in 1:length(Tf)  ] 
end


function is_relevant(mons1, mons2, incidence1, incidence2, dict1, dict2, cells, verts)
    bools = []

    for j in 1:length(cells)
        cell = cells[j]

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
                w1 = (v2-v1)//gcd(v2-v1)
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


function new_ray(Tf, cell, mons1, mons2, incidence1, incidence2, dict1, dict2, two_cells, verts)
    two_cell_indices = filter(index -> issubset(cell, two_cells[index]), 1:length(two_cells))
    j = two_cell_indices[1]
    mon1 = mons1[dict1[j]]
    mon2 = mons2[dict2[j]]

    v1 = reverse(Int64.(mon1.exps))
    v2 = reverse(Int64.(mon2.exps))
    
    inds = findall(i -> verts[cell[i],1] == 0, 1:length(cell))
    res = verts[cell[inds], 2:end] * [v1 v2]
    if rank(res) > 1
        error("Applying Tf to cell gives something positive dimensional.")
    end
    #if rank(res) == 0
    #    error("Applying Tf to cell gives back zero.")
    #end
    if rank(res) == 0
        return zeros(Int64, 2)
    end
    index_nonzero_vector = findfirst(i -> sum(res[i,:].!=0)!=0, 1:size(res, 1))
    ray = res[index_nonzero_vector,:]
    ray = Vector{Int64}(ray/abs(gcd(ray)))
end



function get_rays(V1, V2)
    Tf1, Tf2, Tf = trop_pols(V1, V2)
    TH = TropicalHypersurface(Tf)
    #visualize(TH.polyhedralComplex)

    TH, codim_1_cells, two_cells, verts = top_cells(Tf)
    mons1 = monomials(Tf1)
    incidence1, dict1 = mon_incidence(mons1, two_cells, verts)
    mons2 = monomials(Tf2)
    incidence2, dict2 = mon_incidence(mons2, two_cells, verts)

    is_rel = is_relevant(mons1, mons2, incidence1, incidence2, dict1, dict2, two_cells, verts)
    is_lat = is_lateral(mons1, mons2, incidence1, incidence2, dict1, dict2, two_cells, verts, is_rel)
    is_ess, is_adj = is_essential(mons1, mons2, incidence1, incidence2, dict1, dict2, two_cells, verts, is_lat, is_rel, codim_1_cells)
    #relev_inds = findall(is_rel)

    ray_list = [[0; 1], [1; 0]]

    for i in 1:length(codim_1_cells)
        line  = codim_1_cells[i]
        if is_ess[i] 
            println(i)
            ray = new_ray(Tf, line, mons1, mons2, incidence1, incidence2, dict1, dict2, two_cells, verts)
            push!(ray_list, ray)
        end
    end

    for i in 1:length(codim_1_cells)
        line  = codim_1_cells[i]
        if is_adj[i] 
            ray = new_ray(Tf, line, mons1, mons2, incidence1, incidence2, dict1, dict2, two_cells, verts)
            push!(ray_list, ray)
        end
    end

    for i in 1:length(two_cells)
        cell = two_cells[i]
        if is_rel[i]
            ray = new_ray(Tf, cell, mons1, mons2, incidence1, incidence2, dict1, dict2, two_cells, verts)
            push!(ray_list, ray)
        end
    end
    unique!(ray_list)
    sort!(ray_list, by = ray -> atan(ray[2], ray[1]))
end


function get_ray_mults(V1, V2, ray_list, lpts1, lpts2)    
    @polyvar X Y U V 
    f1 = randn(Float64, length(lpts1))'*[X^pt[1]*Y^pt[2] for pt in lpts1]
    f2 = randn(Float64, length(lpts2))'*[X^pt[1]*Y^pt[2] for pt in lpts2]
    Jac = differentiate(f1, X) * differentiate(f2, Y) - differentiate(f1, Y) * differentiate(f2, X)
    Newt_Jac = hcat(Jac.x.Z ...)'

    dual_ray_list = []
    for ray in ray_list
        push!(dual_ray_list, [-ray[2], ray[1]])
        push!(dual_ray_list, [ray[2], -ray[1]])
    end
    sort!(dual_ray_list, by = ray -> atan(ray[2], ray[1]))
    vecs = [dual_ray_list[i]+dual_ray_list[i+1] for i in 1:length(dual_ray_list)÷2]
    #vecs = [rand(-100:100, 2) for i in 1:round(length(ray_list))+1]
    vecs = [vec.÷gcd(vec) for vec in vecs]

    M_inters = zeros(Int64, length(ray_list), length(vecs))
    is_pos = zeros(Bool, length(ray_list), length(vecs))
    m_vols = zeros(Int64, length(vecs))

    for i in 1:length(ray_list)
        for j in 1:length(vecs)
            vec = vecs[j]
            ray = ray_list[i]
            is_pos[i, j] = ray'*vec > 0
            M_inters[i, j] = abs( ray[2] * vec[2] + ray[1] * vec[1] )
            #abs(det([[-ray[2];ray[1]] vec]))
        end
    end

    for i in 1:length(vecs)
        vec = vecs[i]
        a = [-min(0, vec[1]); -min(0, vec[2])]
        b = vec+a
        pol_C = rand(Float64)*U^a[1]*V^a[2] + rand(Float64)*U^b[1]*V^b[2]
        vol = mixed_volume([f1 - U, f2 - V, Jac, pol_C])
        m_vols[i] = vol
    end

    M = transpose([M_inters.*is_pos M_inters.*(.!is_pos)])
    M = MatrixSpace(ZZ, size(M)...)(M)
    m = [m_vols; m_vols]
    m = MatrixSpace(ZZ, length(m),1)(m)
    if rank(M) < length(ray_list)
        error("The incidence matrix is not of full rank!")
    end
    sol = Oscar.solve(M, m)
    Matrix{Int64}(sol)[:]
end


function reconstr_discr(verts1, verts2, lpts1, lpts2)
    ray_list = get_rays(verts1, verts2)
    ray_mults = get_ray_mults(verts1, verts2, ray_list, lpts1, lpts2) 
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

