


verts1 = transpose([0 2; 2 2; 4 4; 2 6]);
verts2 = transpose([1 2; 2 2; 5 5; 3 6]);

data1 = get_data_fast(verts1)[2];
data2 = get_data_fast(verts2)[2];
Jac = newt_Jac(data1, data2)
face_list = list_faces(data1, data2)
ray_list = get_rays_slim(data1, data2,Jac)
newt_D = cpt_discriminant_chow(data1, data2, Jac, face_list, ray_list)
Delta = get_D_elim(data1.verts, data2.verts)
Bool(Polymake.polytope.congruent(newt_D.pm_polytope, Delta.pm_polytope))





correction = correction_term(data1, data2, face_list)
weight_list = []
test_pol = get_newton_polytope(ray_list)
A = transpose(hcat(ray_list...))
b = get_inequalities(ray_list, test_pol)

test_volume = mixed_volume_via_pullback(data1, data2, Jac, test_pol, correction)
@assert test_volume == volume(Delta + test_pol) - volume(Delta) - volume(test_pol)
for i in 1:length(ray_list)
    b_mod = copy(b)
    b_mod[i] += 1
    other_test_pol = Polyhedron(-A,-b_mod)
    #mixed_vol(correct_newt_D, test_pol) - mixed_vol(correct_newt_D, other_test_pol)

    println("i = $i")
    m_vol = mixed_volume_via_pullback(data1, data2, Jac, other_test_pol, correction) 
    m_vol_correct = volume(Delta + other_test_pol) - volume(Delta) - volume(other_test_pol)

    println("correct_vol = $m_vol")
    println("computed_vol = $m_vol_correct")

    @assert m_vol == m_vol_correct

    weight = test_volume - m_vol
    push!(weight_list, weight)
end









pol = other_test_pol




vol = volume(Jac + pol_lifted) - volume(Jac) - volume(pol_lifted)
(a_max, b_max) = Int64.([ maximum([lpt[1] for lpt in lpts]) , maximum([lpt[2] for lpt in lpts]) ])
vol = (vol - dot(correction, [a_max; b_max]))#[1]

vol = volume(pol1 + pol2) - volume(pol1) - volume(pol2)
a_max = maximum(exps[1,:])
b_max = maximum(exps[2,:])
vol = (vol - dot(correction, [a_max; b_max]))#[1]





verts = vertices(pol)
transl = [ -minimum([vert[1] for vert in verts]) , -minimum([vert[2] for vert in verts]) ]
pol_transl = pol+transl


S,(z1,z2) = PolynomialRing(ZZ,2)

lpts = lattice_points(pol)
transl = [ -minimum([lpts[1] for lpts in lpts]) , -minimum([lpts[2] for lpts in lpts]) ]
lpts = [lpt + transl for lpt in lpts]
h = rand(1:1, length(lpts))'*[z1^pt[1]*z2^pt[2] for pt in lpts]
(f1, f2) = (data1.f, data2.f)
g = h(f1, f2)
pol2 = newton_polytope(g)


pol1 = Jac
pol2_new = lifted_polytope(data1, data2, pol_transl)


vol = volume(pol1 + pol2) - volume(pol1) - volume(pol2)

vol_new = volume(pol1 + pol2_new) - volume(pol1) - volume(pol2_new)
vol = vol_new
if !(vol == vol_new)
    verts_pol = vertices(pol)
    println("verts = $verts_pol")
end
@assert vol == vol_new


exps = hcat(collect(exponents(h))...)
a_max = maximum(exps[1,:])
b_max = maximum(exps[2,:])
vol = (vol - dot(correction, [a_max; b_max]))#[1]

vol2 = Polymake.polytope.mixed_volume(correct_newt_D.pm_polytope, newton_polytope(h).pm_polytope)

if !(vol == QQ(vol2))
    verts_pol = vertices(pol_transl)
    println("verts = $verts_pol")
end
@assert(vol == QQ(vol2))
