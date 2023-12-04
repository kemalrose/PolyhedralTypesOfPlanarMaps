


needsPackage "tropicalToric"

A = QQ[z1,z2,y1,y2]


verts1 = transpose matrix {{1,0},{2,0},{0,3}};
allverts1 = fold( (v, w) -> v|w,  latticePoints convexHull verts1)
f1 = sum apply(entries transpose allverts1, vert -> (random(500)*z1^(vert_0)*z2^(vert_1)));

verts2 = transpose matrix {{1,0},{0,1},{2,0},{1,1}};
allverts2 = fold( (v, w) -> v|w,  latticePoints convexHull verts2)
f2 = sum apply(entries transpose allverts2, vert -> (random(500)*z1^(vert_0)*z2^(vert_1)));


Jac = jacobian(ideal(f1,f2));
Jac = submatrix'(Jac,{2,3},{});
Jideal = minors(2,Jac);

JJac = jacobian(ideal(Jideal_0,f1));
JJac = submatrix'(JJac,{2,3},{});
JJideal = minors(2,JJac);


I = ideal(Jideal, f1-y1, f2-y2);

Dideal = eliminate({z1,z2}, I);
Cideal = eliminate({y1,y2}, I);

(vertsJ,vertsJJ,vertsD,vertsC) = apply((Jideal,JJideal,Dideal,Cideal), ideal -> vertices newtonPolytope ideal_0)
vertsJ = submatrix'(vertsJ, {2,3},{});
vertsJJ = submatrix'(vertsJJ, {2,3},{});
vertsD = submatrix'(vertsD, {0,1},{});
vertsC = submatrix'(vertsC, {2,3},{});

(J,JJ,D,C) = apply((vertsJ,vertsJJ,vertsD,vertsC),verts -> convexHull verts)

A1 = convexHull verts1
A2 = convexHull verts2
A01 = convexHull (verts1 | matrix {{0},{0}})
A02 = convexHull (verts2 | matrix {{0},{0}})

m1 = min(  (entries (submatrix'(verts1,{1},{})|submatrix'(verts2,{1},{})))_0) - 1
m2 = min(  (entries (submatrix'(verts1,{0},{})|submatrix'(verts2,{0},{})))_0) - 1


newtonsNumber = 1
sigma1 = 1
sigma2 = 0

"In the torus:"
lambda1 = degree ideal(  diff(z1,f1) , diff(z2,f1) )
lambda2 = degree ideal(  diff(z1,f2) , diff(z2,f2) )

Deltagamma = 

p1 = mixedVolume {A01,A02}
p2 = length interiorLatticePoints J
p3 = m1+m2 
p4 = m1*m2
p5 = mixedVolume{J,JJ} - lambda1 + sigma2
p6 = (length interiorLatticePoints D) - length interiorLatticePoints J - mixedVolume{J,JJ} + lambda1 + sigma1

(p1, p2, p3, p4, p5, p6)