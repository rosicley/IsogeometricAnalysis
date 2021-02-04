p0 = newp; Point(p0) = {10.000000, 0.125000, 0.0, 1.000000}; Physical Point('p0') = {p0};
//
p1 = newp; Point(p1) = {10.000000, 0.000000, 0.0, 1.000000}; Physical Point('p1') = {p1};
//
p2 = newp; Point(p2) = {20.000000, 0.000000, 0.0, 1.000000}; Physical Point('p2') = {p2};
//
p3 = newp; Point(p3) = {20.000000, 0.125000, 0.0, 1.000000}; Physical Point('p3') = {p3};
//
l0 = newl; Line(l0) = {p0, p1}; Physical Line('l0') = {l0};
//
l1 = newl; Line(l1) = {p1, p2}; Physical Line('l1') = {l1};
//
l2 = newl; Line(l2) = {p2, p3}; Physical Line('l2') = {l2};
//
l3 = newl; Line(l3) = {p3, p0}; Physical Line('l3') = {l3};
//
ll0 = newll; Line Loop(ll0) = {l0, l1, l2, l3};
//
s0 = news; Plane Surface(s0) = {ll0}; Physical Surface('s0') = {s0};
//
Mesh.ElementOrder = 3;
//
Mesh.Algorithm = 2;
//
Mesh 2;
//
Mesh.MshFileVersion = 2.2;
Save "initial.msh";
