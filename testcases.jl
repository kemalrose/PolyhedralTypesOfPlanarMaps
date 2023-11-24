verts1 = [0 1 3 3 1; 3 3 1 0 1]
verts2 = [0 0 1 3 2 1; 1 2 3 1 0 0]
D = get_discriminant(verts1, verts2)
data1, data2, cbd, newt_D, ray_list = get_stuff(verts1, verts2)
#get_ray_mults(cbd, ray_list, data1, data2)



Observe that the monomial terms of the polynomial det Jac f are expressed as (c1 + ... + cr)*x^(a - 1), where
 1:= (1,1) , and each c_i is of the form m*n*det(u,v), 
m*z^u and n*z^u are monomial terms of the respective (generically chosen) polynomials f1 and f2 such that u + v = a .

lpts = lattice_points(convex_hull(verts1)+convex_hull(verts2))
filter(lpt -> , lpts)
1- Take the set S := (A1 + A2) \ K, where

K is the set of all points a in the Minkowski sum A1 + A2
for which any u in A1 and v in A2 with a = u + v satisfies det(u,v) = 0. 

2- Translate S by the vector (-1, -1), to obtain S':= S + {(-1, -1)}

3- Remove all points from S' with a negative coordinate to obtain
S'' := S' \ {(b1, b2) : b1 < 0 or b2 < 0}

4- The Newton polytope P is the convex hull of S''.





verts1 = [0  1  3  3  1;
          3  3  1  0  1]
verts2 = [0  0  1  3  2  1;
          1  2  3  1  0  0]
D = get_discriminant(verts1, verts2)
data1, data2, cbd, newt_D, ray_list = get_stuff(verts1, verts2)
get_ray_mults(cbd, ray_list, data1, data2)
D_newt = get_D_elim(verts1, verts2)


verts1 = 
[0  1  3;
1  2  0]
verts2 = 
[0  1  2  2;
2  2  1  0]
D = get_discriminant(verts1, verts2)
data1, data2, cbd, newt_D, ray_list = get_stuff(verts1, verts2)

D_newt = get_D_elim(verts1, verts2)
#cbd = combine_data(data1, data2, Vars)
#rays_list = get_rays(cbd)

potential_rays = []
for (i,cell) in enumerate(cbd.codim_1_cells )
    n_ray = new_ray(cell, cbd::combined_data)
    n_2_cells =  length(cbd.two_cells)
    two_cell_indices = filter(index -> issubset(cell, cbd.two_cells[index]), 1:n_2_cells)
    #two_cell_1 = cbd.two_cells[two_cell_indices[1]]

    cbd.mons1[cbd.dict1[1]]
    cbd.mons1[cbd.dict1[2]]
    cbd.mons2[cbd.dict2[1]]
    cbd.mons2[cbd.dict2[2]]

    if (n_ray[1] == -3) && (n_ray[2] == -2)
        print(i,cell)
    end
    push!(potential_rays, n_ray)
end
unique!(potential_rays)
hcat(potential_rays...)


verts1 = [  0  0  3;
            1  3  0 ]
verts2 = [0  0  3  2;
1  2  0  0
]
D = get_discriminant(verts1, verts2)

verts1 = [0  0  2; 1  3  1];
verts2 = [0  0  1  2; 1  2  2  1];
D = get_discriminant(verts1, verts2)


verts1 = transpose([0 2; 2 2; 4 4; 2 6]);
verts2 = transpose([1 2; 2 2; 5 5; 3 6]);
D = get_discriminant(verts1, verts2)
psi = get_topological_type(verts1, verts2)


verts1 = hcat([[0,1], [0,3], [3,6], [1,2]]...)
verts2 = hcat([[1,0], [2,0], [1,1], [4,4]]...)
D = get_discriminant(verts1, verts2)
psi = get_topological_type(verts1, verts2)


A1 = (0,1), (0,3), (3,6), (1,2)

A2 = (1,0), (2,0), (1,1), (4,4)

\Psi_1 = 9

\Psi_2 = 9

\Psi_3 = \Psi_4 = 0

\Psi_5 = [I don't have the software to compute this, maybe you can check it using elimination also]

\Psi_6 = -7

\Psi_8 = 12


pol_list = listAllPolytopes(5)
for pol in pol_list
    polyt = convex_hull(pol')
    if (size(pol,2) != length(vertices(polyt)))
        error("wrong representation")
    end
end

polyt = convex_hull(verts2')
lpts = Vector{Vector{Int64}}(lattice_points(polyt))
is_conic(lpts)

verts1 =
[0  1  3;
1  2  0]
verts2 =
[0  1  2  2;
2  2  1  0]