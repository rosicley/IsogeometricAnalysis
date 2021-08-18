p0 = newp; Point(p0) = {0.000000, 0.000000, 0.0, 0.150000}; Physical Point('p0') = {p0};
//
p1 = newp; Point(p1) = {2.000000, 0.000000, 0.0, 0.150000}; Physical Point('p1') = {p1};
//
p2 = newp; Point(p2) = {2.000000, 1.000000, 0.0, 0.150000}; Physical Point('p2') = {p2};
//
p3 = newp; Point(p3) = {2.000000, 2.000000, 0.0, 0.150000}; Physical Point('p3') = {p3};
//
p4 = newp; Point(p4) = {0.000000, 2.000000, 0.0, 0.150000}; Physical Point('p4') = {p4};
//
p5 = newp; Point(p5) = {0.500000, 0.000000, 0.0, 0.100000}; Physical Point('p5') = {p5};
//
p6 = newp; Point(p6) = {0.500000, 0.250000, 0.0, 0.000500}; Physical Point('p6') = {p6};
//
p7 = newp; Point(p7) = {0.500000, 1.500000, 0.0, 0.150000};
//
p8 = newp; Point(p8) = {0.650000, 1.500000, 0.0, 0.150000}; Physical Point('p8') = {p8};
//
p9 = newp; Point(p9) = {0.500000, 1.650000, 0.0, 0.150000}; Physical Point('p9') = {p9};
//
p10 = newp; Point(p10) = {0.350000, 1.500000, 0.0, 0.150000}; Physical Point('p10') = {p10};
//
p11 = newp; Point(p11) = {0.500000, 1.350000, 0.0, 0.150000}; Physical Point('p11') = {p11};
//
p12 = newp; Point(p12) = {1.000000, 1.000000, 0.0, 0.150000};
//
p13 = newp; Point(p13) = {1.000000, 1.150000, 0.0, 0.150000}; Physical Point('p13') = {p13};
//
p14 = newp; Point(p14) = {0.850000, 1.000000, 0.0, 0.150000}; Physical Point('p14') = {p14};
//
p15 = newp; Point(p15) = {1.150000, 1.000000, 0.0, 0.150000}; Physical Point('p15') = {p15};
//
p16 = newp; Point(p16) = {1.000000, 0.850000, 0.0, 0.150000}; Physical Point('p16') = {p16};
//
p17 = newp; Point(p17) = {1.500000, 0.500000, 0.0, 0.150000};
//
p18 = newp; Point(p18) = {1.500000, 0.650000, 0.0, 0.150000}; Physical Point('p18') = {p18};
//
p19 = newp; Point(p19) = {1.350000, 0.500000, 0.0, 0.150000}; Physical Point('p19') = {p19};
//
p20 = newp; Point(p20) = {1.650000, 0.500000, 0.0, 0.150000}; Physical Point('p20') = {p20};
//
p21 = newp; Point(p21) = {1.500000, 0.350000, 0.0, 0.150000}; Physical Point('p21') = {p21};
//
l0 = newl; Line(l0) = {p0, p5}; Physical Line('l0') = {l0};
//
l1 = newl; Line(l1) = {p5, p1}; Physical Line('l1') = {l1};
//
l2 = newl; Line(l2) = {p1, p2}; Physical Line('l2') = {l2};
//
l3 = newl; Line(l3) = {p2, p3}; Physical Line('l3') = {l3};
//
l4 = newl; Line(l4) = {p3, p4}; Physical Line('l4') = {l4};
//
l5 = newl; Line(l5) = {p4, p0}; Physical Line('l5') = {l5};
//
l6 = newl; Circle(l6) = {p8, p7, p9}; Physical Line('l6') = {l6};
//
l7 = newl; Circle(l7) = {p9, p7, p10}; Physical Line('l7') = {l7};
//
l8 = newl; Circle(l8) = {p10, p7, p11}; Physical Line('l8') = {l8};
//
l9 = newl; Circle(l9) = {p11, p7, p8}; Physical Line('l9') = {l9};
//
l10 = newl; Circle(l10) = {p15, p12, p13}; Physical Line('l10') = {l10};
//
l11 = newl; Circle(l11) = {p13, p12, p14}; Physical Line('l11') = {l11};
//
l12 = newl; Circle(l12) = {p14, p12, p16}; Physical Line('l12') = {l12};
//
l13 = newl; Circle(l13) = {p16, p12, p15}; Physical Line('l13') = {l13};
//
l14 = newl; Circle(l14) = {p20, p17, p18}; Physical Line('l14') = {l14};
//
l15 = newl; Circle(l15) = {p18, p17, p19}; Physical Line('l15') = {l15};
//
l16 = newl; Circle(l16) = {p19, p17, p21}; Physical Line('l16') = {l16};
//
l17 = newl; Circle(l17) = {p21, p17, p20}; Physical Line('l17') = {l17};
//
p22 = newp; Point(p22) = {0.500000, 0.215000, 0.0, 0.100000};
//
p23 = newp; Point(p23) = {0.535000, 0.250000, 0.0, 0.100000};
//
p24 = newp; Point(p24) = {0.500000, 0.285000, 0.0, 0.100000};
//
p25 = newp; Point(p25) = {0.465000, 0.250000, 0.0, 0.100000};
//
c0 = newl;
 //
c00 = newl; Line(c00) = {p5, p22};
//
c01 = newl; Line(c01) = {p22, p6};
//
Physical Line('c0') = {c00, c01}; 
//
j00 = newl; Circle(j00) = {p22, p6, p23}; 
//
j01 = newl; Circle(j01) = {p23, p6, p24}; 
//
j02 = newl; Circle(j02) = {p24, p6, p25}; 
//
j03 = newl; Circle(j03) = {p25, p6, p22}; 
//
Physical Line('j0') = {j00, j01, j02, j03};
//
ll0 = newll; Line Loop(ll0) = {l0, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17};
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
Plugin(Crack).OpenBoundaryPhysicalGroup = p5;
Plugin(Crack).Run; 
//
Mesh.MshFileVersion = 2.2;
Save "initial.msh";
