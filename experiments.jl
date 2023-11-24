
using Oscar, HomotopyContinuation, Primes, LinearAlgebra, Random

include("listPolygons.jl")
include("cpt_D.jl")
include("computePsiUsingOscar.jl")


A1, A2 = [5 5 0; 1 0 0], [6 1 0; 0 5 0];
show(15*A1)
show(15*A2)
A1, A2 = 
@time Delta = get_delta(A1, A2);

A1, A2 = [15 15 0; 3 0 0], [18 3 0; 0 15 0];
@time Delta = get_delta(A1, A2);



show(Primes.prime.(1:100))
#f1 = 75*z1^5*z2 + 15*z1^5 + 95*z1^2 + 35*z1 + 258
#f2 = 250*z1^6 + 136*z1^5*z2 + 253*z1^5 + 63*z1^4*z2^2 +   160*z1*z2^5 +43*z1*z2^2 + 168*z1*z2 + 64*z1 + 63
    p = 139
    FF,_ = FiniteField(p,1,"")
    R,(z1,z2,y1,y2) = PolynomialRing(FF,["z1","z2","y1","y2"]);
    f1 = 41*z1^5*z2 + 43*z1^5 + 47*z1^2 + 53*z1 + 59
    f2 = 61*z1^6 + 67*z1^5*z2 + 71*z1^5 + 73*z1^4*z2^2 +   79*z1*z2^5 + 83*z1*z2^2 + 89*z1*z2 + 97*z1 + 101   
    Jac = det(jacobi_matrix([f1,f2])[1:2,:]);
    I = ideal(R,[Jac, f1-y1, f2-y2]);
    @time Discriminant = eliminate(I,[z1,z2])[1];
    D = convex_hull([vert[3:4] for vert in vertices(newton_polytope(D_elim))])


    vol = volume(D)
    println("p = $p")
    println("vol = $vol")
    println(hcat(vertices(D)...)) 
    verts1 = hcat(vertices(newton_polytope(f1))...)[1:2,:]
    verts2 = hcat(vertices(newton_polytope(f2))...)[1:2,:]
    verts1 = Int64.(verts1)
    verts2 = Int64.(verts2)

    verts1 = [5 5 0; 1 0 0]
    verts2 = [6 1 0; 0 5 0]
    data1 = get_data_fast(verts1)[2];
    data2 = get_data_fast(verts2)[2];
    
    Jac = newt_Jac(data1, data2)
    face_list = list_faces(data1, data2)
    ray_list = get_rays_slim(data1, data2,Jac)
    newt_D = cpt_discriminant_chow(data1, data2, Jac, face_list, ray_list)
    volume(newt_D)
    
    A1, A2 = [5 5 0; 1 0 0], [6 1 0; 0 5 0];
    @time Delta = get_delta(A1, A2);




    
    p = Primes.prime(rand(1:100))
        #p = 541 is an interesting choice
        #just kidding
    FF,_ = FiniteField(p,1,"")
    R,(z1,z2,y1,y2) = PolynomialRing(FF,["z1","z2","y1","y2"]);
    f1 = 41*z1^5*z2 + 43*z1^5 + 47*z1^2 + 53*z1 + 59
    f2 = 61*z1^6 + 67*z1^5*z2 + 71*z1^5 + 73*z1^4*z2^2 +   79*z1*z2^5 + 83*z1*z2^2 + 89*z1*z2 + 97*z1 + 101    


    #f1 = 75*z1^5*z2 + 15*z1^5 + 95*z1^2 + 35*z1 + 258
    #f2 = 250*z1^6 + 136*z1^5*z2 + 253*z1^5 + 63*z1^4*z2^2 +   160*z1*z2^5 +43*z1*z2^2 + 168*z1*z2 + 64*z1 + 63    

    verts1 = hcat(vertices(newton_polytope(f1))...)[1:2,:]
    verts2 = hcat(vertices(newton_polytope(f2))...)[1:2,:]
    #random_poly( convex_hull(transpose(verts1)) ,R);
    #f1 = random_poly( convex_hull(transpose(verts1)) ,R);
    #f2 = random_poly(convex_hull(transpose(verts2)),R);

    #random_poly( convex_hull(verts1) ,R);
    #f2 = random_poly(convex_hull(verts2'),R);
    Jac = det(jacobi_matrix([f1,f2])[1:2,:]);
    I = ideal(R,[Jac, f1-y1, f2-y2]);
    D_elim = eliminate(I,[z1,z2])[1];
    D = convex_hull([vert[3:4] for vert in vertices(newton_polytope(D_elim))])
    vol = volume(D)
    println("p = $p")
    println("vol = $vol")
    println(hcat(vertices(D)...))

    QQFieldElem[49 19 0 0; 0 30 30 0]


    f1 = random_poly( convex_hull(transpose(verts1)) ,R);
    f2 = random_poly(convex_hull(transpose(verts2)),R);
    Jac = det(jacobi_matrix([f1,f2])[1:2,:]);
    I = ideal(R,[Jac, f1-y1, f2-y2]);
    D_elim = eliminate(I,[z1,z2])[1];
    D = convex_hull([vert[3:4] for vert in vertices(newton_polytope(D_elim))])
    vol = volume(D)
    println("p = $p")
    println("vol = $vol")
    println(hcat(vertices(D)...))

    69  45  16  0   0
    0   40  50  50  0
   

    lattice_points(convex_hull(transpose(verts1)))
    lattice_points(convex_hull(transpose(verts2)))
    
    visualize(convex_hull(transpose(verts1)))
    visualize(convex_hull(transpose(verts2)))


    f1 = 75*z1^5*z2 + 15*z1^5 + 95*z1^2 + 35*z1 + 258
f2 = 250*z1^6 + 136*z1^5*z2 + 253*z1^5 + 63*z1^4*z2^2 +   160*z1*z2^5 + 176*z1*z2^4 + 195*z1*z2^3 + 43*z1*z2^2 + 168*z1*z2 + 64*z1 + 63


f1 = 75*z1^5*z2 + 15*z1^5 + 95*z1^2 + 35*z1 + 258
f2 = 250*z1^6 + 136*z1^5*z2 + 253*z1^5 + 63*z1^4*z2^2 +   160*z1*z2^5 +43*z1*z2^2 + 168*z1*z2 + 64*z1 + 63



verts1 = Int64.(verts1)
verts2 = Int64.(verts2)
data1 = get_data_fast(verts1)[2];
data2 = get_data_fast(verts2)[2];
Jac = newt_Jac(data1, data2)
face_list = list_faces(data1, data2)
ray_list = get_rays_slim(data1, data2,Jac)
newt_D = cpt_discriminant_chow(data1, data2, Jac, face_list, ray_list)

hcat(vertices(newt_D)...)
volume(newt_D)

Delta = get_D_elim(data1.verts, data2.verts)
@assert Polymake.polytope.congruent(newt_D.pm_polytope, Delta.pm_polytope) == 1









pol_Lists =  list_all_polytopes.([7,8])

pol1 = convex_hull([[0,0],[0,0],[0,0],[0,0]])
pol2 = convex_hull([[0,0],[0,0],[0,0],[0,0]])



Primes.prime(5)
Primes.prime(7)
Primes.prime(100)
Primes.prime(108)
Primes.prime(57)
Primes.prime(135)

for i in 1:10
    p = Primes.prime(rand(1:100))
        #p = 541 is an interesting choice

    FF,_ = FiniteField(p,1,"")
    R,(z1,z2,y1,y2) = PolynomialRing(FF,["z1","z2","y1","y2"]);
    f1 =  11 + 541*z1^10*z2 +  1223 * z2^4
    f2 =  761 + 17*z1^6 +  593 * z1^4*z2^5 
    #random_poly( convex_hull(verts1') ,R);
    #f2 = random_poly(convex_hull(verts2'),R);
    Jac = det(jacobi_matrix([f1,f2])[1:2,:]);
    I = ideal(R,[Jac, f1-y1, f2-y2]);
    D_elim = eliminate(I,[z1,z2])[1];
    D = convex_hull([vert[3:4] for vert in vertices(newton_polytope(D_elim))])
    vol = volume(D)
    println("p = $p")
    println("vol = $vol")
end


random_poly( convex_hull(verts1') ,R)
f1 = random_poly( convex_hull(verts1') ,R);
f2 = random_poly(convex_hull(verts2'),R);




verts1 = [2 2 4 3; 1 2 0 0]
verts2 = [0 0 1 2 1; 2 3 3 1 1]



verts1 = hcat([[5, 1], [0, 4], [0, 0]]...)
verts2 = hcat([[6, 0], [4, 3], [0, 0]]...)
data1 = get_data_fast(verts1)[2];
data2 = get_data_fast(verts2)[2];
Jac = newt_Jac(data1, data2)
face_list = list_faces(data1, data2)
ray_list = get_rays_slim(data1, data2,Jac)
newt_D = cpt_discriminant_chow(data1, data2, Jac, face_list, ray_list)
Delta = get_D_elim(data1.verts, data2.verts)
@assert Polymake.polytope.congruent(newt_D.pm_polytope, Delta.pm_polytope) == 1











function do_certification(deg::Int)
   
    type_list = []
    data_list = list_all_conical_pols(deg)


    for index in 1:200
    #for i in 1:length(data_list)
    #    for j in i:length(data_list)
                i = rand(1:length(data_list))
                j = rand(1:length(data_list))
                (data1,data2) = data_list[[i,j]]
                Jac = newt_Jac(data1, data2)
                face_list = list_faces(data1, data2)
                ray_list = get_rays_slim(data1, data2,Jac)


                #correct_newt_D = get_D_elim(data1.verts, data2.verts);

                newt_D = cpt_discriminant_chow(data1, data2, Jac, face_list, ray_list)
                correct_newt_D = get_D_elim(data1.verts, data2.verts);
                @assert Polymake.polytope.congruent(newt_D.pm_polytope, correct_newt_D.pm_polytope) == 1
                
                
                #get_D(data1, data2, Jac, face_list)

                #D = cpt_D_from_scratch(data1.verts, data2.verts);
                #bool = !(sort(rays(normal_fan(D))) == sort(rays(normal_fan(D_newt))))
                if !Bool(Polymake.polytope.congruent(newt_D.pm_polytope, correct_newt_D.pm_polytope))
                    correct_newt_D = get_D_elim(data1.verts, data2.verts);
                    if !Bool(Polymake.polytope.congruent(newt_D.pm_polytope, correct_newt_D.pm_polytope))
                        println("Get error for:")
                        println("verts1 = $verts1")
                        println("verts2 = $verts2")

                        correct_rays = [r.//gcd(r) for r in rays(normal_fan(correct_newt_D))]
                        correct_rays = [Int64[ray[i].num for i in 1:length(ray)] for ray in correct_rays]
                        rays_mults_correct = (ray_list ==  sort!(correct_rays, by = ray -> atan(ray[2], ray[1])))

                        correct_mults = Int64[]
                        Lpts_correct = Oscar.lattice_points(correct_newt_D)
                        for ray in correct_rays
                            min_val = minimum([(transpose(ray)*l_pt) for l_pt in Lpts_correct])
                            mult = length(findall(l_pt -> (transpose(ray)*l_pt) == min_val, Lpts_correct) )
                            push!(correct_mults, mult-1)
                        end
                        rays_correct = (ray_mults == correct_mults)
                        println("The rays are correct: $rays_correct")
                        println("The ray mults are correct: $rays_mults_correct")

                        println("Computed ray mults : $ray_mults")
                        println("Correct ray mults : $correct_mults")

                        error("")
                    end
                else
                    #println("Asserted correctness for index =", index )
                    println("Asserted correctness for index =", (i,j) )
                end

        #end
    end
    type_list
end




n_pols = length.(list_all_polytopes.([1,2,3,4,5,6,7,8]))
length.(list_all_conical_pairs([1,2,3,4,5]))
data_list = list_all_conical_pairs(5)
types = do_experiments(3)

con_list = list_all_conical_pols(5);


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












A1 = [0 1 0 3;
      1 0 3 0]
A2 = A1
verts1 = A1
verts2 = A2
data1 = get_data_fast(verts1)[2]
data2 = get_data_fast(verts2)[2]
visualize(data1.pol)
visualize(data2.pol)
Δ = get_D(A1, A2)
ψ = get_polyhedral_type(A1,A2)

Delta = get_D_elim(data1.verts, data2.verts)


Jac = newt_Jac(data1, data2)
face_list = list_faces(data1, data2)
ray_list = get_rays_slim(data1, data2,Jac)
D = cpt_discriminant_chow(data1, data2, Jac, face_list, ray_list)
correct_newt_D = get_D_elim(data1.verts, data2.verts)
mv, mh = get_m(data1,data2)
N0 = get_Newton_nr(D);


face_list = list_faces(data1, data2)

is_dicritical.(face_list)









A1 = 2*[1 2 3 2; 1 2 0 0]
A2 = 2*[0 0 1 4 3; 2 3 3 0 0]
A1, A2 = verts1, verts2

A1 = [0 2 4 2; 2 2 4 6];
A2 = [1 2 5 3; 2 2 5 6];
Δ = get_D(A1, A2)
ψ = get_polyhedral_type(A1,A2)

verts1 = [1 2 3 2; 1 2 0 0]
verts2 = [0 0 1 4 3; 2 3 3 0 0]


verts1 = [1 2 3 2; 1 2 0 0]
verts2 = [0 0 1 4 3; 2 3 3 0 0]



verts1 = [2 2 4 3; 1 2 0 0]
verts2 = [0 0 1 2 1; 2 3 3 1 1]
verts1 = [0 0 3 1; 1 3 1 0]
verts2 = [0 0 1 2; 2 3 3 2]

verts1 = [2 2 4 3; 1 2 0 0]
verts2 = [0 0 1 2 1; 2 3 3 1 1]


verts1 = [2 2 4 3; 1 2 0 0]
verts2 = [0 0 1 2 1; 2 3 3 1 1]


verts1 = [2 2 4 3; 1 2 0 0]
verts2 = [0 0 1 2 1; 2 3 3 1 1]
The rays are correct: false
The ray mults are correct: false
Computed ray mults : [1, 8, 1, 12, 2, 12]
Correct ray mults : [16, 1, 19, 1, 22]



verts1 = [2 2 4 3; 1 2 0 0]
verts2 = [0 0 1 2 1; 2 3 3 1 1]
The rays are correct: false
The ray mults are correct: false
Computed ray mults : [1, 8, 1, 12, 2, 12]
Correct ray mults : [2, 12, 4, 17, 1, 17]

data1 = get_data_fast(verts1)[2]
data2 = get_data_fast(verts2)[2]
Jac = newt_Jac(data1, data2)
face_list = list_faces(data1, data2)
ray_list = get_rays_slim(data1, data2,Jac)
newt_D = cpt_discriminant_chow(data1, data2, Jac, face_list, ray_list)


ray_mults = get_ray_mults_fast(Jac, ray_list, data1, data2, face_list)

newt_D = pol_from_rays_and_mults(ray_list, ray_mults)
correct_newt_D = get_D_elim(data1.verts, data2.verts);
Polymake.polytope.congruent(convex_hull(transpose(newt_D)).pm_polytope, correct_newt_D.pm_polytope)

correct_verts = sort(vertices(correct_newt_D))

#(sort(vertices(newt_D)) == correct_verts)
correct_rays = [r.//gcd(r) for r in rays(normal_fan(correct_newt_D))]
correct_rays = [Int64[ray[i].num for i in 1:length(ray)] for ray in correct_rays]
ray_list ==  sort!(correct_rays, by = ray -> atan(ray[2], ray[1]))

correct_mults = Int64[]
Lpts_correct = Oscar.lattice_points(correct_newt_D)
for ray in correct_rays
    min_val = minimum([(transpose(ray)*l_pt) for l_pt in Lpts_correct])
    mult = length(findall(l_pt -> (transpose(ray)*l_pt) == min_val, Lpts_correct) )
    push!(correct_mults, mult-1)
end
ray_mults == correct_mults












D  = get_discriminant(verts1, verts2)
@var X Y Z W
f1 = sum([round(randn(),digits = 2)*X^verts1[1,i]*Y^verts1[2,i] for i in 1:size(verts1,2)])
f2 = sum([round(randn(),digits = 2)*X^verts2[1,i]*Y^verts2[2,i] for i in 1:size(verts2,2)])
J = det(HomotopyContinuation.jacobian(System([f1,f2])))
pts = []
for i in 1:250
    l = randn(4)'* [X;Y;Z;W]
    sols = HomotopyContinuation.solve(System([f1-Z,f2-W,J,l]))
    new_pts = [pt[3:4] for pt in solutions(sols)]
    pts = [pts; new_pts]
end
filtered_pts = filter(pt -> abs(norm(pt)) < 1e4, pts)


function get_vandermonde(f1,f2,D)
    @var X Y
    mons = Vector{Int64}.(lattice_points(D))
    N = 1.2*length(mons)
    sample_D(f1, f2, N)
    get_vdm(pts, mons; normalize = true)
end



    

V = get_vdm(filtered_pts, mons)
rank(V) == minimum(size(V))
@time U, sings, VV = svd(V)
        coeff = VV[:, end]
        #println(round.(log10.(sings)))
        err_tol = sings[end]/sings[end-1] * 100
        println(sings[end])
        println(sings[end-1])
        if sings[end] == 0.0
            err_tol = eps(T)/sings[end-1] * 100
        end



@time U, sings, VV = svd(V)
coeff = VV[:, end]
#println(round.(log10.(sings)))
err_tol = sings[end]/sings[end-1] * 100
println(sings[end])
println(sings[end-1])

function get_vdm(pts, mons)
    T = eltype(pts[1])
    V = zeros(T, (length(pts), length(mons)))
    for i in 1:length(pts)
        V[i, :] = [prod(pts[i].^mon) for mon in mons]
        V[i, :] = V[i, :]/norm(V[i, :])
    end
    V
end



function get_samples_D(F1,F2)
    
end


A1, A2 = data1.pol, data2.pol;
f1 = random_poly(A1,R);
f2 = random_poly(A2,R);
Jac = det(jacobi_matrix([f1,f2])[1:2,:]);
I = ideal(R,[Jac, f1-y1, f2-y2]);
C = eliminate(I,[y1,y2])[1];
D = eliminate(I,[z1,z2])[1];

vec = vecs[i]
b = [-min(0, vec[1]); -min(0, vec[2])]
a = vec+b
pol_C = rand(-50:50)*y1^a[1]*y2^a[2] + rand(-50:50)*y1^b[1]*y2^b[2]

A,(z1,z2) = PolynomialRing(FF,["z1","z2"]);
f1 = random_poly(A1,A)
f2 = random_poly(A2,A)
Jac = det(jacobi_matrix([f1,f2]));
g = rand(-50:50)*f1+rand(-50:50)*f2;
J5 = ideal(A,[Jac, g]);
MV = Polymake.polytope.mixed_volume(newton_polytope(Jac).pm_polytope,newton_polytope(g).pm_polytope)
MV - 2*a[1]
MV - a[1]*(l1+l2) - b[2]*(h1+h2)
J5torus = (J5:ideal(A,z1^300)):ideal(A,z2^300);
vdim((quo(A,J5torus)[1]))

#verts_list = list_all_polytopes(4)


cbd = combine_data(data1, data2, Vars)
face_list = list_faces(data1, data2)

D = get_D(data1,data2,cbd,Vars)

hcat(vertices(newton_polytope(f1))...)
hcat(vertices(newton_polytope(f2))...)
hcat(vertices(newton_polytope(D))...)
mons_cbd_Jac = hcat(vertices(convex_hull([mon.z for mon in cbd.Jac.x]))...)
hcat(vertices(newton_polytope(det(jacobi_matrix([f1,f2])[1:2,:])))...)

#data_list = []
#for verts in verts_list
#    V = transpose(verts)
#    is_con, Data = get_auxillary_data(V,verts,RR,Vars)
#    if is_con
#        push!(data_list, Data)
#    end
#end
#
#data1, data2 = data_list[[end-1, end]]
#@time cbd = combine_data(data1, data2, Vars);
#D_verts = get_D_verts(cbd)



p = 1299709
p = Primes.prime(1000200)
FF,_ = FiniteField(p,1,"")
R,(Z1,Z2,Y1,Y2) = PolynomialRing(FF,["Z1","Z2","Y1","Y2"]);

A1, A2 = data1.pol, data2.pol;
F1 = random_poly(A1,R);
F2 = random_poly(A2,R);
Jac = det(jacobi_matrix([F1,F2])[1:2,:]);
I = ideal(R,[Jac, F1-Y1, F2-Y2]);
C = eliminate(I,[Y1,Y2])[1];
DD = eliminate(I,[Z1,Z2])[1];
VD = vertices(newton_polytope(DD))
newt_D = convex_hull(transpose(hcat(VD...)[3:4,:]))

binomial = rand(-50:50)*Y1^a[2]+rand(-50:50)*Y2^b[1]
G = ideal(R,[binomial, Jac, F1-Y1, F2-Y2])
Gtorus = G:ideal(R,(Z1*Z2)^100)
vdim(quo(R,Gtorus)[1])
H = ideal(R,[DD,binomial, Z1, Z2])
mixed_volume

length(interior_lattice_points(newt_D))

@time cbd = combine_data(data1, data2);
@time ray_list = get_rays(cbd);
@time ray_mults = get_ray_mults(cbd, ray_list, data1, data2);    
@time res = pol_from_rays_and_mults(ray_list, ray_mults);




V1 = vcat([[0 2], [2 2], [4 4], [2 6]]...);
V2 = vcat([[1 2], [2 2], [5 5], [3 6]]...);

V1 = [0 2; 2 0]
V2 = [1 0; 2 2; 3 1]
V_res = [0 2; 12 2; 0 8; 1 0; 10 0]
verts1 = transpose(V1)		
verts2 = transpose(V2)
bool1, data1 = get_auxillary_data(V1, verts1, RR, Vars)
bool2, data2 = get_auxillary_data(V2, verts1, RR, Vars)
lpts1, lpts2 = data1.lpts, data2.lpts
cbd = combine_data(data1, data2)
ray_list = get_rays(cbd)
ray_mults = get_ray_mults(cbd, ray_list, data1, data2)    
res = pol_from_rays_and_mults(ray_list, ray_mults)


verts_corr = vertices(newton_polytope(D))
P_res = Oscar.convex_hull([vert[3:4] for vert in verts_corr])
correct_rays = [r.//gcd(r) for r in rays(normal_fan(P_res))]
correct_rays = [Int64[ray[i].num for i in 1:length(ray)] for ray in correct_rays]
sort!(correct_rays, by = ray -> atan(ray[2], ray[1]))

correct_mults = Int64[]
Lpts_correct = Oscar.lattice_points(P_res)
for ray in correct_rays
    min_val = minimum([(transpose(ray)*l_pt) for l_pt in Lpts_correct])
    mult = length(findall(l_pt -> (transpose(ray)*l_pt) == min_val, Lpts_correct) )
    push!(correct_mults, mult-1)
end




V1 = [0 1; 0 2; 2 1; 2 2]
V2 = [0 2; 2 1; 1 2; 2 0]		
V_res = [0 1; 11 0; 8 4; 2 0; 0 12]
verts1, verts2 = transpose(V1), transpose(V2) 
ray_list = get_rays(verts1, verts2)
facets(convex_hull(V_res))

lpts1, lpts2 = [Vector{Vector{Int64}}(lattice_points(convex_hull(verts))) for verts in [verts1, verts2]]
ray_mults = get_ray_mults(verts1, verts2, ray_list, lpts1, lpts2) 







A1, A2 = convex_hull(verts1), convex_hull(verts2);

ψ1, ψ2, ψ3, ψ4, ψ5, ψ6, ψ7, ψ8 = combinatorial_type(A1, A2)

p = 1299709
FF,_ = FiniteField(p,1,"")
A,(z1,z2,y1,y2) = PolynomialRing(FF,["z1","z2","y1","y2"])
RR,(Z1,Z2) = PolynomialRing(QQ,["Z1","Z2"])
f1 = random_poly(A1,R);
f2 = random_poly(A2,R);

allverts1 = Vector{Int64}.(lattice_points(convex_hull(verts1)))
coeffsf1 = rand(-15:15, length(allverts1))'
f1 = coeffsf1*[z1^pt[1]*z2^pt[2] for pt in allverts1]
F1 = coeffsf1*[Z1^pt[1]*Z2^pt[2] for pt in allverts1]

allverts2 = Vector{Int64}.(lattice_points(convex_hull(verts2)))
coeffsf2 = rand(-15:15, length(allverts2))'
f2 = coeffsf2*[z1^pt[1]*z2^pt[2] for pt in allverts2]
F2 = coeffsf2*[Z1^pt[1]*Z2^pt[2] for pt in allverts2]

Jac = det(jacobi_matrix([f1,f2])[1:2,:])
Jideal = ideal(A,Jac)

JJac = det(jacobi_matrix([Jac,f1])[1:2,:]);
JJideal = ideal(A,JJac);

I = ideal(A,[Jac, f1-y1, f2-y2]);
Cideal = eliminate(I,[y1,y2]);
Dideal = eliminate(I,[z1,z2]);
D = Dideal[1]

ideal_D_sings = ideal(A,[derivative(D,y1),derivative(D,y2),D,z1,z2]):ideal(A,[y1^120,y2^120]);
coord_ring = quo(A,ideal_D_sings)[1];
nr_sings_D = vdim(coord_ring)





( npolJ,npolJJ,npolD,npolC) = [newton_polytope(ideal[1]) for ideal in [Jideal, JJideal, Dideal, Cideal]]
(vertsJ,vertsJJ,vertsD,vertsC) = [transpose(hcat(vertices(pol)...)) for pol in [npolJ,npolJJ,npolD,npolC]]
vertsJ = vertsJ[:,1:2]
vertsJJ = vertsJJ[:,1:2]
vertsD = vertsD[:,3:4]
vertsC = vertsC[:,1:2]
(A1,A2,J,JJ,D,C) = convex_hull.([verts1,verts2,vertsJ,vertsJJ,vertsD,vertsC])


id1 = (ideal(A,[derivative(f1,z1),derivative(f1,z2),y1,y2]):ideal(A,[z1])):ideal(A,[z2])
id2 = (ideal(A,[derivative(f2,z1),derivative(f2,z2),y1,y2]):ideal(A,[z1])):ideal(A,[z2])
AI1,_ = quo(A, id1); AI2,_ = quo(A, id2);
λ1 = vdim(AI1); λ2 = vdim(AI2)

m1 = minimum( hcat([verts1;verts2]...)[1,:] )
m2 = minimum( hcat([verts1;verts2]...)[2,:] )


visualize(A1);
visualize(A2);
visualize(J);
visualize(JJ);
visualize(D);
visualize(C);
Σ = normal_fan(convex_hull(verts1)+ convex_hull(verts2))
visualize(Σ)
face_list = list_faces(verts1, verts2)

for face in face_list
    println("")
    (γ1,γ2,ω) = face
    println("γ1 = $γ1")
    println("γ2 = $γ2")
    println("ω = $ω")
    bool_dicritical = is_dicritical(γ1,γ2,ω)
    bool_long = is_long(γ1,γ2,ω)
    bool_origin = is_origin(γ1,γ2,ω)
    bool_semi_origin = is_semi_origin(γ1,γ2,ω)
    bool_lower = is_lower(γ1,γ2,ω)
    bool_upper = is_upper(γ1,γ2,ω)
    println("dicritical = $bool_dicritical ")
    println("long = $bool_long ")
    println("origin = $bool_origin ")
    println("semi_origin = $bool_semi_origin ")
    println("lower = $bool_lower ")
    println("upper = $bool_upper ")
end



function LA(verts1, verts2, allverts1, allverts2, coeffsf1, coeffsf2, face_list)
    if sum([is_dicritical(γ1, γ2, ω) for (γ1, γ2, ω) in face_list]) == 0
        0
    end
    lower_long =[]
    for (γ1, γ2, ω) in face_list
        if (is_lower(γ1, γ2, ω) && is_long(γ1, γ2, ω))
            push!(lower_long, (γ1, γ2, ω))
        end
    end

    lower_long_not_origin =[]
    for (γ1, γ2, ω) in lower_long
        if !is_origin(γ1, γ2, ω)
            push!(lower_long_not_origin, (γ1, γ2, ω))
        end
    end

    if length(lower_long) == 0
        return -1
    end
    if length(lower_long_not_origin) > 0
        (γ1, γ2, ω) = lower_long_not_origin[1]
    else
        (γ1, γ2, ω) = lower_long[1]
    end
    f1ω, min1 = restrictfct([allverts1; [[0; 0]]], [coeffsf1 [-y1]], ω)
    f2ω, min2 = restrictfct([allverts2; [[0; 0]]], [coeffsf2 [-y2]], ω)

    IGω = ideal(A,[f1ω,f2ω])
    IGω = saturation(saturation(IGω,ideal(A,z1)),ideal(A,z2))
    IGω = eliminate(IGω,[z1,z2])$$›
    verts = [vert[3:4] for vert in vertices(newton_polytope(IGω[1]))]
    Δγ = convex_hull(verts)
    nlptsΔγ =  length(interior_lattice_points(Δγ))
    if length(lower_long_not0) > 0
        return nlptsΔγ  + 2
    else
        return nlptsΔγ - σ(1,Δγ) + 1
end


verts = [vert[3:4] for vert in vertices(newton_polytope(IGω[1]))]
Δγ = convex_hull(verts)

            


p1 = mixed_vol( [verts1;[[0,0]]], [verts2;[[0,0]]] )
p2 = length(interior_lattice_points(J))
p3 = m1+m2 
p4 = m1*m2
p5 = mixed_vol(vertsJ, vertsJJ) - λ1 + σ2(D)
p6 = length(interior_lattice_points(D)) - length(interior_lattice_points(J)) - mixed_vol(vertsJ, vertsJJ) + λ1 + σ1(D)
p7 = 
p8 = 




f_0 = 37\,x\,y+21\,x\,z+22\,y\,z+47\,x+2\,y+7\,z+35,\\
f_1 = 13\,x^{5}y^{5}z^{5}+49\,x+44\,y+2\,z+14


f1 = 4\,z_1^{4}\,z_2^{4} + 8\,z_1^{3}\,z_2^{5} + 9\,z_1^{3}\,z_2^{4} - z_1^{3}\,z_2^{3} + 15\,z_1^{2}\,z_2^{6} - 15\,z_1^{2}\,z_2^{5} - 14\,z_1^{2}\,z_2^{4} - 13\,z_1^{2}\,z_2^{3} - 11\,z_1^{2}\,z_2^{2} - 12\,z_1\,z_2^{4} - 4\,z_1\,z_2^{3} - 10\,z_1\,z_2^{2} - 14\,z_2^{2}
f2 = 10\,z_1^{5}\,z_2^{5} + 5\,z_1^{4}\,z_2^{5} - 5\,z_1^{4}\,z_2^{4} + 7\,z_1^{3}\,z_2^{6} - 5\,z_1^{3}\,z_2^{5} + 14\,z_1^{3}\,z_2^{4} - 12\,z_1^{3}\,z_2^{3} - 4\,z_1^{2}\,z_2^{4} - 15\,z_1^{2}\,z_2^{3} + 14\,z_1^{2}\,z_2^{2} - 15\,z_1\,z_2^{2}



