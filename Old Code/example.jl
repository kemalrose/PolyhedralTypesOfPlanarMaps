using Oscar



A,(x,y) = PolynomialRing(QQ,["x","y"]);
f = x^2-y^3
f = x * y
sing = ideal(A,[f,derivative(f,x),derivative(f,y)])
vdim(quo(A,sing:ideal(A,[x,y]))[1])


function random_poly(newton_poly, R)
    allverts = Vector{Int64}.(lattice_points(newton_poly))
    coeffs = rand(-50:50, length(allverts))'
    coeffs*[R[1]^pt[1]*R[2]^pt[2] for pt in allverts]
end



p = 1299709
p = Primes.prime(1000100)
FF,_ = FiniteField(p,1,"")
R,(z1,z2,y1,y2) = PolynomialRing(FF,["z1","z2","y1","y2"]);
verts1 = [[0, 2], [2, 2], [4, 4], [2, 6]];
verts2 = [[1, 2], [2, 2], [5, 5], [3, 6]];

verts1 = [[0, 2], [2, 2], [3, 3], [1, 4]];
verts2 = [[1, 2], [2, 2], [3, 3], [2, 4]];

A1, A2 = convex_hull(verts1), convex_hull(verts2);
f1 = random_poly(A1,R);
f2 = random_poly(A2,R);

fiber_ideal = ideal(R, [f1 - rand(-50:50),f2 - rand(-50:50),y1,y2]);
degreeF = vdim(quo(R,fiber_ideal)[1])

Jac = det(jacobi_matrix([f1,f2])[1:2,:]);
I = ideal(R,[Jac, f1-y1, f2-y2]);
C = eliminate(I,[y1,y2])[1];
D = eliminate(I,[z1,z2])[1];

interior_lattice_points(newton_polytope(Jac))
interior_lattice_points(newton_polytope(D))


ψ5 = 
ψ6 = 
2*ψ5 + ψ6 + mult_of_0




Rnew,(z1,z2,t) = GradedPolynomialRing(QQ,["z1","z2","t"]);
coeffs = transpose([c.n for c in collect(coefficients(C))])
Cstar = coeffs * [z1^exp[1] * z2^(exp[2]-3) * t^(16-exp[1]-exp[2]) for exp in exponents(C)]
Cstart = Rnew(Cstar)
Crve = Oscar.ProjPlaneCurve(Cstar)
is_smooth_curve(Crve)
arithmetic_genus(Crve)
geometric_genus(Crve)
length(interior_lattice_points(newton_polytope(Jac)))
sing_locus_c_star = ideal(Rnew, [Cstar,derivative(Cstar,z1),derivative(Cstar, z2),t])
minimal_generators(sing_locus_c_star)
dim(sing_locus_c_star)
sing_pol = eliminate(sing_locus_c_star, [z1])


Dstar = transpose(collect(coefficients(D))) * [z1^exp[1] * z2^(exp[2]) * t^(36-exp[1]-exp[2]) for exp in exponents(D)]
CrveD = ProjPlaneCurve(Dstar)
is_smooth_curve(CrveD)
arithmetic_genus(CrveD)
geometric_genus(CrveD)






vdim(quo(R,sing_locus_c_star)[1])


D = eliminate(I,[z1,z2])[1];
ideal_D_sings = ideal(R,[derivative(D,y1),derivative(D,y2),D,z1,z2]);
ideal_D_sings_at_zero = ideal(R,[derivative(D,y1),derivative(D,y2),D,z1,z2,y1^100,y2^100]);
radical_ideal_D_sings = Oscar.radical(ideal_D_sings);

sings_D = radical(ideal_D_sings)
#sings_D has 5 gens
pol_y = eliminate(sings_D, [y1])
#Nr of singularities
nr_sings_w_o_mults = vdim(quo(R,radical_ideal_D_sings)[1]) - 1
nr_sings_w_mults = vdim(quo(R,ideal_D_sings)[1])
nr_sings_w_mults_at_zero = vdim(quo(R,ideal_D_sings_at_zero)[1])

interior_lattice_points(newton_polytope(Jac))
interior_lattice_points(newton_polytope(D))
interior_lattice_points(newton_polytope(C))




degD = 35


#Arithmetic genus: 
g_arithmetic  = 1//2*(degD-1)*(degD-2)

#Geometrix genus:
g_geometric =  nr_sings - g_arithmetic

det(jacobi_matrix([D])[1:2,:]);
ideal(R, )