
CFLAGS	         = -I${PETSC_DIR}/src/dm/mesh/sieve ${BOOST_INCLUDE} ${TRIANGLE_INCLUDE} ${TETGEN_INCLUDE}
FFLAGS	         = 
CPPFLAGS         = 
FPPFLAGS         =
LOCDIR           = src/ksp/ksp/examples/tutorials/
MANSEC           = KSP
CLEANFILES       = rhs.vtk solution.vtk
NP               = 1
DISPLAY          = :0.0

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

petsc_poisson: petsc_poisson.o  chkopts
	-${CLINKER} -o petsc_poisson petsc_poisson.o  ${PETSC_KSP_LIB}
	${RM} petsc_poisson.o

petsc_poisson_v2: petsc_poisson_v2.o chkopts
	 -${CLINKER} -o petsc_poisson_v2 petsc_poisson_v2.o ${PETSC_KSP_LIB}
	${RM} petsc_poisson_v2.o

petsc_poisson_par: petsc_poisson_par.o chkopts
	-${CLINKER} -o petsc_poisson_par petsc_poisson_par.o ${PETSC_KSP_LIB}
	${RM} petsc_poisson_par.o

petsc_poisson_par_barry: petsc_poisson_par_barry.o chkopts
	-${CLINKER} -o petsc_poisson_par_barry petsc_poisson_par_barry.o ${PETSC_KSP_LIB}
	${RM} petsc_poisson_par_barry.o

petsc_poisson_par_barry2: petsc_poisson_par_barry2.o chkopts
	-${CLINKER} -o petsc_poisson_par_barry2 petsc_poisson_par_barry2.o ${PETSC_KSP_LIB}
	${RM} petsc_poisson_par_barry2.o



#----------------------------------------------------------------------------
run_petsc_poisson:
	-@${MPIEXEC} -n 1 ./petsc_poisson 


include ${PETSC_DIR}/conf/test
