MANSEC           = KSP
CLEANFILES       = rhs.vtk solution.vtk
NP               = 1
CSOURCES 				 = $(wildcard *.cpp Solid/BoundaryConditions/*.cpp Solid/Isogeometric/*.cpp Solid/*.cpp Solid/FiniteElement/*.cpp Solid/Mesh/*.cpp)
FCOMPILER        = gfortran -O2
CXXFLAGS        += -w
CURRENT_DIR      = $(shell pwd)

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

s:$(GCH_FILES:.h=.gch) $(CSOURCES:.cpp=.o)
	@-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB} -lboost_system -std=c++0x 

debug: $(CSOURCES:.cpp=.o)
	@-${CLINKER} -o $@ $^ ${PETSC_KSP_LIB} -g
	@gdb debug

clear:
	@$ rm *.o *~ s *.vtu mirror* ISOdomain* FEdomain* *.mod *.dat ma26* tensao* esforc* saida omega.txt *.geo Solid/BoundaryConditions/*.o Solid/Isogeometric/*.o Solid/*.o Solid/FiniteElement/*.o Solid/Mesh/*.o $(CURRENT_DIR)/src/*.gch

run1:
	@$ mpirun -np 1 ./s  

run2:
	@$ mpirun -np 2 ./s

run3:
	@$ mpirun -np 3 ./s

run4teste:
	@$ mpirun -n 4 ./s -ksp_view -ksp_monitor_singular_value

run5teste:
	@$ mpirun -n 5 ./s -ksp_view -pc_type none -ksp_type gmres -ksp_monitor_singular_value -ksp_gmres_restart 1000

run4:
	@$ mpirun -n 4 ./s

run5:
	@$ mpirun -np 5 ./s

run6:
	@$ mpirun -np 6 ./s 

run7:
	@$ mpirun -np 7 ./s

run8:
	@$ mpirun -np 8 ./s

run12:
	@$ mpirun -np 12 ./s

run10:
	@$ mpirun -np 10 ./s

run11:
	@$ mpirun -np 11 ./s

run6teste:
	@$ mpirun -np 6 ./s -pc_type cholesky -pc_factor_mat_solver_type mumps -ksp_view -mat_mumps_icntl_4 2 -mat_mumps_icntl_5 0 -mat_mumps_icntl_18 3 -mat_mumps_icntl_28 2 -mat_mumps_icntl_29 2 -ksp_type gmres -ksp_monitor_singular_value -ksp_gmres_restart 1000

run6teste2:
	@$ mpirun -np 6 ./s -ksp_view -pc_type cholesky -ksp_type gmres -ksp_monitor_singular_value -ksp_gmres_restart 1000
