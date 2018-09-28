static char help[] = "Reads a PETSc matrix and vector from a socket connection,  solves a linear system and sends the result back.\n";

/*T
   Concepts: KSP^solving a linear system
   Processors: n
T*/

/* 
  Include "petscksp.h" so that we can use KSP solvers.  Note that this file
  automatically includes:
     petscsys.h    - base PETSc routines   petscvec.h - vectors
     petscmat.h - matrices
     petscis.h     - index sets        petscksp.h - Krylov subspace methods
     petscviewer.h - viewers           petscpc.h  - preconditioners
*/
#include "petscksp.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  KSP            ksp;          /* linear solver context */
  PC             prec;         // preconditioner 
  //MatNullSpace   nsp;
  Mat            A;            /* matrix */
  Vec            x,b,test;     /* approx solution, RHS, exact solution */
  PetscViewer    fd;           /* viewer */
  PetscErrorCode ierr;
  PetscScalar    val;
  PetscReal      rnorm;
  PetscInt       its;
  PetscLogDouble time1,time2;

  PetscInitialize(&argc,&args,(char *)0,help);
  fd = PETSC_VIEWER_SOCKET_WORLD;
  ierr = VecCreateSeq(PETSC_COMM_SELF,1,&test);CHKERRQ(ierr);

  ierr = MatLoad(fd,MATAIJ,&A);CHKERRQ(ierr);
  
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER);CHKERRQ(ierr);
  // CG
  ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);
  ierr = KSPCGSetType(ksp,KSP_CG_SYMMETRIC);
  // ICC
  KSPGetPC(ksp,&prec);
  PCSetType(prec,PCICC);
  KSPSetTolerances(ksp,1.e-14,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSetUp(ksp);CHKERRQ(ierr);
  
  while (0<1)
  {
      ierr = VecLoad(fd,VECSEQ,&b);CHKERRQ(ierr);

      ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
      ierr = VecSet(x,0);CHKERRQ(ierr);

      ierr = PetscGetTime(&time1);CHKERRQ(ierr);

      ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
      ierr = KSPGetIterationNumber(ksp,&its);
      ierr = KSPGetResidualNorm(ksp,&rnorm); 

      ierr = PetscGetTime(&time2);CHKERRQ(ierr);
      val = time2-time1;


      // send to matlab
      ierr = VecView(x,fd);CHKERRQ(ierr);
      ierr = VecSet(test,its);CHKERRQ(ierr);
      ierr = VecView(test,fd);CHKERRQ(ierr);
      ierr = VecSet(test,rnorm);CHKERRQ(ierr); 
      ierr = VecView(test,fd);CHKERRQ(ierr);
      //ierr = VecSet(test,val);CHKERRQ(ierr); 
      //ierr = VecView(test,fd);CHKERRQ(ierr); 
  
  }    
  // destroy vectors and matrices
  ierr = MatDestroy(A);CHKERRQ(ierr); 
  ierr = VecDestroy(b);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = KSPDestroy(ksp);CHKERRQ(ierr); 

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

