p0 = newp; Point(p0) = {0.000000, 0.000000, 0.0, 1.000000}; Physical Point('p0') = {p0};
//
p1 = newp; Point(p1) = {10.000000, 0.000000, 0.0, 1.000000}; Physical Point('p1') = {p1};
//
p2 = newp; Point(p2) = {10.000000, 10.000000, 0.0, 1.000000}; Physical Point('p2') = {p2};
//
p3 = newp; Point(p3) = {0.000000, 10.000000, 0.0, 1.000000}; Physical Point('p3') = {p3};
//
p4 = newp; Point(p4) = {4.566987, 4.750000, 0.0, 0.050000}; Physical Point('p4') = {p4};
//
p5 = newp; Point(p5) = {5.433013, 5.250000, 0.0, 0.050000}; Physical Point('p5') = {p5};
//
p6 = newp; Point(p6) = {0.000000, 5.000000, 0.0, 1.000000}; Physical Point('p6') = {p6};
//
l0 = newl; Line(l0) = {p0, p1}; Physical Line('l0') = {l0};
//
l1 = newl; Line(l1) = {p1, p2}; Physical Line('l1') = {l1};
//
l2 = newl; Line(l2) = {p2, p3}; Physical Line('l2') = {l2};
//
l3 = newl; Line(l3) = {p3, p6}; Physical Line('l3') = {l3};
//
l4 = newl; Line(l4) = {p6, p0}; Physical Line('l4') = {l4};
//
p7 = newp; Point(p7) = {5.216506, 5.125000, 0.0, 0.050000}; Physical Point('p7') = {p7};
//
p8 = newp; Point(p8) = {5.558013, 5.033494, 0.0, 0.050000}; Physical Point('p8') = {p8};
//
p9 = newp; Point(p9) = {5.649519, 5.375000, 0.0, 0.050000}; Physical Point('p9') = {p9};
//
p10 = newp; Point(p10) = {5.308013, 5.466506, 0.0, 0.050000}; Physical Point('p10') = {p10};
//
c0 = newl;
 //
c00 = newl; Line(c00) = {p4, p7};
//
c01 = newl; Line(c01) = {p7, p5};
//
Physical Line('c0') = {c00, c01}; 
//
j00 = newl; Circle(j00) = {p7, p5, p8}; 
//
j01 = newl; Circle(j01) = {p8, p5, p9}; 
//
j02 = newl; Circle(j02) = {p9, p5, p10}; 
//
j03 = newl; Circle(j03) = {p10, p5, p7}; 
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
Mesh.ElementOrder = 2;
//
Mesh.Algorithm = 2;
//
Mesh 2;
//
Plugin(Crack).Dimension = 1;
Plugin(Crack).PhysicalGroup = c0;
Plugin(Crack).Run; 
//
Mesh.MshFileVersion = 2.2;
Save "initial.msh";
