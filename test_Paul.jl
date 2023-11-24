
  
  verts1 = transpose([0 2; 2 2; 4 4; 2 6]);
  verts2 = transpose([1 2; 2 2; 5 5; 3 6]);
  verts1 = [0  1  3  3  1;
  3  3  1  0  1]
verts2 = [0  0  1  3  2  1;
  1  2  3  1  0  0]
  verts1 = [0 1 3 3 1; 3 3 1 0 1]
  verts1 = [0 0 1 3 2 1; 1 2 3 1 0 0]


@time get_rays_slim(verts1, verts2)

"""
  RR = TropicalSemiring(min)
    S,(x,y) = RR["x","y"]
    p = Primes.prime(1000200)
    FF,_ = FiniteField(p,1,"")
    R,(z1,z2) = PolynomialRing(FF,["z1","z2"]);
    Vars = (x,y,z1,z2)

V1 = transpose(verts1)
V2 = transpose(verts2)
data1 = get_auxillary_data(V1, verts1, RR, Vars)[2];
data2 = get_auxillary_data(V2, verts2, RR, Vars)[2];
cbd = combine_data(data1,data2,Vars)
ray_list = get_rays(cbd)
"""

