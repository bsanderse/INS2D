static char help[] = "Reads a PETSc matrix and vector from a socket connection,  solves a linear system and sends the result back.\n";

/*T
   Concepts: KSP^solving a linear system
   Processors: n
T*/

/* 
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h       - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
*/
#include "petscksp.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{

  Vec            b;        /* approx solution, RHS, exact solution */
  PetscViewer    fd;              /* viewer */
  PetscErrorCode ierr;


  PetscInitialize(&argc,&args,(char *)0,help);
  
  fd   = PETSC_VIEWER_SOCKET_WORLD;

 
  // load rhs vector
  ierr = VecLoad(fd,VECMPI,&b);CHKERRQ(ierr);

    // send to matlab
    ierr = VecView(b,fd);CHKERRQ(ierr);
    ierr = VecDestroy(b);CHKERRQ(ierr);
    PetscFinalize();CHKERRQ(ierr);
  return 0;
}

