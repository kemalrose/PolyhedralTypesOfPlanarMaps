





X = normal_toric_variety(correct_newt_D)
CRing = Oscar.chow_ring(X)
Id = CRing.I
X_rays = rays(X)

class = sum(gens(CRing).*b)
dot(gens(CRing),gens(CRing))


coeffs = []
verts = [Array(vert) for vert in vertices(correct_newt_D)]
for ray in X_rays
    ray = ray .* lcm(denominator.(ray))
    weight = -minimum([dot(ray, vert) for vert in verts])
    push!(coeffs, weight)
end
class = sum([coeffs[i]*A[i] for i in 1:length(X_rays)])

class*class
volume(correct_newt_D)*2






verts1 = hcat([[1,1], [3,3], [2,3], [1,2]]...)
verts2 = hcat([[1,1], [2,1], [3,2], [3,3]]...)

verts1 = [1 2 3 2; 1 2 0 0]
verts2 = [0 0 1 4 3; 2 3 3 0 0]

data1 = get_data_fast(verts1)[2]
data2 = get_data_fast(verts2)[2]

new_psi = topological_type(data1,data2)


Jac = newt_Jac(data1, data2)
face_list = list_faces(data1, data2)
ray_list = get_rays_slim(data1, data2,Jac)
#ray_mults = get_ray_mults_fast(Jac, ray_list, data1, data2, face_list)



newt_D = get_discriminant(verts1, verts2)
correct_newt_D = get_D_elim(data1.verts, data2.verts);
@assert Polymake.polytope.congruent(newt_D.pm_polytope, correct_newt_D.pm_polytope) == 1

