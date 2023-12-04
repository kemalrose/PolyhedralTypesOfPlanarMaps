

for p in Primes.primes(1,300)
  FF,_ = FiniteField(p,1,"");
  R,(z1,z2,y1,y2) = PolynomialRing(FF,["z1","z2","y1","y2"]);
  f1 = 41*z1^5*z2 + 43*z1^5 + 47*z1^2 + 53*z1 + 59;
  f2 = 61*z1^6 + 67*z1^5*z2 + 71*z1^5 + 73*z1^4*z2^2
  + 79*z1*z2^5 + 83*z1*z2^2 + 89*z1*z2 + 97*z1+ 101;
  Jac = det(jacobi_matrix([f1,f2])[1:2,:]);
  Id = Oscar.ideal(R,[Jac, f1-y1, f2-y2]);
  Discriminant = eliminate(Id,[z1,z2])[1];
  D = convex_hull([vert[3:4] for vert in vertices(newton_polytope(Discriminant))])
  v = volume(D)
  println("For p = $p the volume is equal to $v")
end


A1, A2 = [5 5 0; 1 0 0], [6 1 0; 0 5 0];
verts1 = [[5,1], [5,0], [0,0]]
verts2 = [[6,0], [1,5], [0,0]]
for p in Primes.primes(300,400)
  FF,_ = FiniteField(p,1,"");
  R,(z1,z2,y1,y2) = PolynomialRing(FF,["z1","z2","y1","y2"]);
  f1 = 108*z1^5*z2 + 66*z1^5 + 64*z1^4 + 85*z1^3 + 53*z1 + 132;
  f2 = 61*z1^6 + 67*z1^5*z2 + 71*z1^5 + 73*z1^4*z2^2
  + 79*z1*z2^5 + 83*z1*z2^2 + 89*z1*z2 + 97*z1+ 101;
  Jac = det(jacobi_matrix([f1,f2])[1:2,:]);
  Id = Oscar.ideal(R,[Jac, f1-y1, f2-y2]);
  Discriminant = eliminate(Id,[z1,z2])[1];
  D = convex_hull([vert[3:4] for vert in vertices(newton_polytope(Discriminant))])
  v = volume(D)
  println("For p = $p the volume is equal to $v")




  D_fin  = compute_Discriminant_over_finite_field(A1, A2, FF)
end







A1 = Int64.(hcat([vert[1:2] for vert in vertices(newton_polytope(f1))]...))
A2 = Int64.(hcat([vert[1:2] for vert in vertices(newton_polytope(f2))]...))
@time Delta = get_delta(A1, A2);
volume(Delta)

vertices(newton_polytope(f2))
A1, A2 = [5 5 0; 1 0 0], [6 1 0; 0 5 0];
@time Delta = get_delta(A1, A2);
volume(Delta)






vertices(newton_polytope(f1))
vertices(newton_polytope(f2))
A1, A2 = [5 5 0; 1 0 0], [6 1 0; 0 5 0];
@time Delta = get_delta(A1, A2);
volume(Delta)


pol_Lists =  list_all_polygons.([1,2,3,4,5])
conic_pols = [filter(verts->is_conic(verts), pol_list) for pol_list in pol_Lists]
length.(conic_pols)  



A1 = [0 2 4 2; 2 2 4 6];
A2 = [1 2 5 3; 2 2 5 6];
Delta = get_delta(A1, A2);
psi = get_polyhedral_type(A1,A2);



  verts1 = transpose([0 2; 2 2; 4 4; 2 6]);
  verts2 = transpose([1 2; 2 2; 5 5; 3 6]);
  data1 = get_auxillary_data(verts1)
  data2 = get_auxillary_data(verts2)


  verts1 = [0  1  3  3  1;
  3  3  1  0  1]
verts2 = [0  0  1  3  2  1;
  1  2  3  1  0  0]
  verts1 = [0 1 3 3 1; 3 3 1 0 1]
  verts1 = [0 0 1 3 2 1; 1 2 3 1 0 0]


@time get_rays_slim(verts1, verts2)


test_nr = 30
deg = 4
pol_list = list_all_polygons(deg)
conic_pols = filter(verts->is_conic(verts), pol_list)

#data_list = get_auxillary_data.(conic_pols)

for h in 1:test_nr
    index1 = rand(1:length(conic_pols))
    index2 = rand(1:length(conic_pols))
    (verts1, verts2) = conic_pols[[index1, index2]]
    Discriminant = get_delta(verts1,verts2)
    Discriminantp = compute_discriminant_probailistically(verts1, verts2)
    if volume(Discriminant)!=volume(Discriminantp)
      println("Something went wrong:")
      println("verts1 = $verts1")
      println("verts2 = $verts2")
    end
end


for h in 1:test_nr
  index1 = rand(1:length(conic_pols))
  index2 = rand(1:length(conic_pols))

  (verts1, verts2) = conic_pols[[index1, index2]]

  data1, data2 =   get_auxillary_data.([verts1, verts2])
  psi = polyhedral_type(data1, data2)
  print(psi)
  if volume(Discriminant)!=volume(Discriminantp)
    println("Something went wrong:")
    println("verts1 = $verts1")
    println("verts2 = $verts2")
  end
end










function compute_discriminant_probailistically(verts1, verts2)
  p1, p2 = Primes.prime(100000), Primes.prime(100000)
  FF1,_ = FiniteField(p1,1,"")
  FF2,_ = FiniteField(p2,1,"")
  Discriminant1 = compute_Discriminant_over_finite_field(verts1, verts2, FF1)
  Discriminant2 = compute_Discriminant_over_finite_field(verts1, verts2, FF2)
  @assert volume(Discriminant1)==volume(Discriminant2)
  Discriminant1
end





function compute_Discriminant_over_finite_field(verts1, verts2, FF)
  R,(z1,z2,y1,y2) = PolynomialRing(FF,["z1","z2","y1","y2"]);

  newton_poly1 = convex_hull(verts1')
  newton_poly2 = convex_hull(verts2')

  allverts1 = Vector{Int64}.(lattice_points(newton_poly1))
  coeffs1 = rand(-100:100, length(allverts1))'
  f1 = coeffs1*[R[1]^pt[1]*R[2]^pt[2] for pt in allverts1]

  allverts2 = Vector{Int64}.(lattice_points(newton_poly2))
  coeffs2 = rand(-100:100, length(allverts2))'
  f2 = coeffs2*[R[1]^pt[1]*R[2]^pt[2] for pt in allverts2]

  Jac = det(jacobi_matrix([f1,f2])[1:2,:]);
  Id = Oscar.ideal(R,[Jac, f1-y1, f2-y2]);
  D_elim = eliminate(Id,[z1,z2])[1];
  D_elim = convex_hull([vert[3:4] for vert in vertices(newton_polytope(D_elim))])
  D_elim
end







