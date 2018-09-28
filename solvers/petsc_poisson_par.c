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

  Mat            A;        /* matrix */
  Vec            b, x;     /* RHS, solution vector */
  KSP            ksp;      /* linear solver context */
  PC             prec;     /* preconditioner */

  PetscViewer    fd;       /* viewer */
  PetscErrorCode ierr;


  PetscInitialize(&argc,&args,(char *)0,help);
  
  // viewer for matlab communication
  fd   = PETSC_VIEWER_SOCKET_WORLD;

  // load matrix
  ierr = MatLoad(fd,MATMPIAIJ,&A);CHKERRQ(ierr);

  // set up krylov space method
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER);CHKERRQ(ierr);
  // CG
  ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);
  ierr = KSPCGSetType(ksp,KSP_CG_SYMMETRIC);CHKERRQ(ierr);
  // ICC
  ierr = KSPGetPC(ksp,&prec);CHKERRQ(ierr);
  ierr = PCSetType(prec,PCBJACOBI);CHKERRQ(ierr);
  //ierr = PCSetType(prec,PCICC);CHKERRQ(ierr);
  //ierr = PCSetType(prec,PCHYPRE);CHKERRQ(ierr);
  //ierr = PCHYPRESetType(prec,"pilut");CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);

  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  
  
  while (0<1)
  {
      // load rhs vector
      ierr = VecLoad(fd,VECMPI,&b);CHKERRQ(ierr);
      
      // duplicate rhs vector to get initialization of solution vector
      ierr = VecDuplicate(b,&x);CHKERRQ(ierr);

      // solve system
      ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

      // send to matlab
      ierr = VecView(x,fd);CHKERRQ(ierr);
  
  }
  
  // destroy vectors and matrices
  ierr = MatDestroy(A);CHKERRQ(ierr); 
  ierr = VecDestroy(b);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = KSPDestroy(ksp);CHKERRQ(ierr); 
  

  PetscFinalize();CHKERRQ(ierr);
  return 0;
}

