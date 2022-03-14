# IsogeometricAnalysis

Isogeometric analysis plays an important role in solid mechanics, offering the possibility of integrating computer-aided design and structural analysis. However, isogeometric discretizations are not as flexible as standard finite elements ones, especially in the context of crack propagation analysis with continuous and complex remeshing. This numerical tool presents a domain decomposition technique for introducing and propagating cracks in an isogeometrically discretized solid. In this technique, a standard triangular finite element mesh is overlapped to the isogeometric model and the basis functions of both discretizations are modified and blended over a region of the physical domain, leading to a new space of functions. This technique allows the insertion and propagation of discontinuities over an isogeometric discretization without modifying the original mesh. The proposed approach is applied to large and small displacements 2D linear elastic fracture mechanics problems, demonstrating to be a very robust, accurate and versatile tool.

##

This project was my [master's work](http://sistemas.set.eesc.usp.br/producao/1250) and resulted in papers published in international conferences and one in [indexed journal](https://www.sciencedirect.com/science/article/pii/S0045782522000366?dgcid=coauthor).
