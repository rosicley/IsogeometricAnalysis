p0 = newp; Point(p0) = {0.000000, 0.000000, 0.0, 0.000200}; Physical Point('p0') = {p0};
//
p1 = newp; Point(p1) = {0.002000, 0.000000, 0.0, 0.000200}; Physical Point('p1') = {p1};
//
p2 = newp; Point(p2) = {0.002000, 0.006000, 0.0, 0.000200}; Physical Point('p2') = {p2};
//
p3 = newp; Point(p3) = {0.000000, 0.006000, 0.0, 0.000200}; Physical Point('p3') = {p3};
//
p4 = newp; Point(p4) = {0.000000, 0.003000, 0.0, 0.000100}; Physical Point('p4') = {p4};
//
p5 = newp; Point(p5) = {0.001000, 0.003000, 0.0, 0.000005}; Physical Point('p5') = {p5};
//
l0 = newl; Line(l0) = {p0, p1}; Physical Line('l0') = {l0};
//
l1 = newl; Line(l1) = {p1, p2}; Physical Line('l1') = {l1};
//
l2 = newl; Line(l2) = {p2, p3}; Physical Line('l2') = {l2};
//
l3 = newl; Line(l3) = {p3, p4}; Physical Line('l3') = {l3};
//
l4 = newl; Line(l4) = {p4, p0}; Physical Line('l4') = {l4};
//
p6 = newp; Point(p6) = {0.000700, 0.003000, 0.0, 0.000050};
//
p7 = newp; Point(p7) = {0.001000, 0.002700, 0.0, 0.000050};
//
p8 = newp; Point(p8) = {0.001300, 0.003000, 0.0, 0.000050};
//
p9 = newp; Point(p9) = {0.001000, 0.003300, 0.0, 0.000050};
//
c0 = newl;
 //
c00 = newl; Line(c00) = {p4, p6};
//
c01 = newl; Line(c01) = {p6, p5};
//
Physical Line('c0') = {c00, c01}; 
//
j00 = newl; Circle(j00) = {p6, p5, p7}; 
//
j01 = newl; Circle(j01) = {p7, p5, p8}; 
//
j02 = newl; Circle(j02) = {p8, p5, p9}; 
//
j03 = newl; Circle(j03) = {p9, p5, p6}; 
//
Physical Line('j0') = {j00, j01, j02, j03};
//
ll0 = newll; Line Loop(ll0) = {l0, l1, l2, l3, l4};
//
s0 = news; Plane Surface(s0) = {ll0}; Physical Surface('s0') = {s0};
//
Line{c00} In Surface{s0};
//
Line{c01} In Surface{s0};
//
Line{j00} In Surface{s0};
//
Line{j01} In Surface{s0};
//
Line{j02} In Surface{s0};
//
Line{j03} In Surface{s0};
//
Mesh.ElementOrder = 3;
//
Mesh.Algorithm = 2;
//
Mesh 2;
//
Plugin(Crack).Dimension = 1;
Plugin(Crack).PhysicalGroup = c0;
Plugin(Crack).OpenBoundaryPhysicalGroup = p4;
Plugin(Crack).Run; 
//
Mesh.MshFileVersion = 2.2;
Save "RJ.msh";
