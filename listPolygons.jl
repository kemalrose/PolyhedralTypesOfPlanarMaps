




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
