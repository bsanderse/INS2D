#include "fintrf.h"
C======================================================================
#if 0
C     
C     cg.F
C     .F file needs to be preprocessed to generate .for equivalent
C     
#endif
C     
C     cg.f
C
C     solve A x = y, with initial guess x_in and tolerance tol
C
C     Fortran:
C     cg(N, A, y, tol, x_in, x, iter, norm1, norm2)
C     Matlab call:
C     [x,iter,norm1,norm2]=cg(A,int64(dia),int64(nd),p,tol,int64(N),x_in,int64(maxit));
C
C     Fortran file with preconditioned conjugate gradient method 
C     linked to Matlab.
C     Original fortran file (7-diagonal matrix): S. v.d. Pijl
C     Extension to general spd matrix: B. Sanderse
C     MEX file: B. Sanderse, April 2010
      
C     This is a MEX-file for MATLAB.
C     Copyright 1984-2009 The MathWorks, Inc.


C======================================================================
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer*8 nlhs, nrhs  

C     Function declarations:
      mwPointer mxCreateNumericArray
      integer*8 mxClassIDFromClassName
      mwPointer mxCreateDoubleMatrix, mxGetPr
      integer*8 mxIsDouble
      mwPointer mxIsInt64
      mwSize mxGetM, mxGetN

C     Pointers to input/output mxArrays:
      mwPointer N, nd, iter
      mwPointer dia, maxiter
      mwPointer A, y, tol, x_in, x_out, norm1, norm2
C     Array information:
      mwSize mrows, ncols, size, one

C     Arguments for mxCreateNumericArray
      integer*8 classid, complexflag
      mwSize ndim, dims(2)

C-----------------------------------------------------------------------
CC    Check for proper number of arguments. 
CC      if(nrhs .ne. 10) then
CC         call mexErrMsgIdAndTxt ('MATLAB:cg:nInput',
CC     +             'Wrong number of input arguments.')
CC      endif
CC    Check for proper number of arguments. 
CC      if(nlhs .ne. 4) then
CC         call mexErrMsgIdAndTxt ('MATLAB:cg:nInput',
CC     +                           'Wrong number of output arguments.')
CC      endif

CC    Validate inputs
      if(mxIsDouble(prhs(1)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:cg:NonNumeric',
     +                           'Input matrix A must be double.')
      endif
      if(mxIsInt64(prhs(2)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:cg:WrongType',
     +     'dia must be int64 vector.')
      endif
      if(mxIsInt64(prhs(3)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:cg:WrongType',
     +     'nd must be int64 scalar.')
      endif
      if(mxIsDouble(prhs(4)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:cg:NonNumeric',
     +                           'Input vector y must be double.')
      endif
      if(mxIsDouble(prhs(5)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:cg:NonNumeric',
     +                           'Input value tol must be double.')
      endif
      if(mxIsInt64(prhs(6)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:cg:WrongType',
     +     'N must be int64 scalar.')
      endif
      if(mxIsDouble(prhs(7)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:cg:NonNumeric',
     +                           'Input vector x must be double.')
      endif
      if(mxIsInt64(prhs(8)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:cg:WrongType',
     +     'maxiter must be int64 scalar.')
      endif
CC    input arguments
C     Create Fortran pointer for input array A
      A = mxGetPr(prhs(1))
C     Create Fortran pointer for input array dia
      dia = mxGetPr(prhs(2))
C     Create Fortran pointer for input array nd
      nd = mxGetPr(prhs(3))
C     Get the size of the input array y
      mrows = mxGetM(prhs(4))
      ncols = mxGetN(prhs(4))
      size = mrows*ncols
C     Create Fortran pointer for input argument y
      y = mxGetPr(prhs(4))
C     Create Fortran pointer for tol.
      tol = mxGetPr(prhs(5))
C     Create Fortran pointer for N.
      N = mxGetPr(prhs(6))
C     Create Fortran pointer for input argument x
      x_in = mxGetPr(prhs(7))
C     Create Fortran pointer for input argument maxiter
      maxiter = mxGetPr(prhs(8))

CC    return arguments
      one = 1
C     Create matrix for the return argument x
      plhs(1) = mxCreateDoubleMatrix(size,one,0)
      x_out = mxGetPr(plhs(1))
C     Create 1x1 integer matrix for the return argument iter
      classid = mxClassIDFromClassName('int64')
      complexflag = 0; ndim = 1; dims(1) = 1
      plhs(2) = mxCreateNumericArray(ndim, dims, 
     &                               classid, complexflag)
      iter = mxGetPr(plhs(2))
C     Create 1x1 matrix for the return argument norm1     
      plhs(3) = mxCreateDoubleMatrix(one,one,0)
      norm1 = mxGetPr(plhs(3))      
C     Create 1x1 matrix for the return argument norm2
      plhs(4) = mxCreateDoubleMatrix(one,one,0)
      norm2 = mxGetPr(plhs(4))

CC    Call the computational subroutine.
      call cg(%val(N), %val(A), %val(dia),  
     &   %val(nd), %val(y), %val(tol), %val(x_in), %val(x_out), 
     &   %val(iter), %val(norm1), %val(norm2), %val(maxiter))

      return
      end

C======================================================================
C     Computational routines

C-------------------------------------------------------------------------------
C  solve for x: Ax = y
C-------------------------------------------------------------------------------
      subroutine cg(N, A, dia, nd, y, tol, x_in, x, 
     &          iter, norm1, norm2, maxiter)
C     general number of diagonals, specified in d
C     integer*8 for input/outputs
C     write replaced by mexErrMsgTxt, stop deleted

      implicit none
      
      integer*8 N, nd, dia(nd)
      double precision A(N,nd), y(N) 
      double precision x_in(N), x(N)      ! input, output     
      integer*8 iter, maxiter
      double precision norm1, norm2, tol

      double precision D(N)
      double precision, dimension(:), allocatable :: b, r_old, z
      double precision, dimension(:), allocatable :: r, p, temp_vec
      
      double precision alpha, temp
      double precision tol2
      double precision rz, rz_old
      double precision norm_res2, norm_res2_0
      
      double precision maxx, minx, maxres
            
      integer*8 i
            
      character*120 line
      integer*8 k, mexPrintf

      allocate(p(N))
      allocate(b(N))
      allocate(r(N))
      allocate(r_old(N))
      allocate(temp_vec(N))
      allocate(z(N))
      
      tol2 = tol*tol

C     if there is only one (main) diagonal the solution can be directly found
      if (nd.eq.1 .and. dia(nd).eq.0) then
        x = y(1:N)/A(1:N,1)
        iter = 1
        norm1 = 0.d0
        norm2 = 0.d0
        goto 110
      end if

C     fill x with x_in
      x  = x_in
     
C     set up preconditioner
      call make_cholinc(N, A, dia, nd, D)
      call matvec_prod(N, A, dia, nd, x, temp_vec)

      r = y - temp_vec

      call precond(N, A, dia, nd, D, r, z)
      
      p = z

      norm_res2_0 = sum(y*y)
      norm_res2   = norm_res2_0

      if (norm_res2_0.le.tol2) then
        iter = 0
	norm1 = 0.d0
	norm2 = 0.d0
	goto 110
      end if
      
      rz = sum(r*z)


      do i=1,maxiter
	
	call matvec_prod(N, A, dia, nd, p, temp_vec)

        temp  = sum(p*temp_vec)
	alpha = rz/temp
	x     = x + alpha*p
	r     = r - alpha*temp_vec
	norm_res2 = sum(r*r)
	norm_res2 = norm_res2 / norm_res2_0
      
        maxx   = maxval(x)
        minx   = minval(x)
        maxres = maxval(abs(r))
	
	if (norm_res2 .le. tol2) goto 100
	
	call precond(N, A, dia, nd, D, r, z)
	
	rz_old = rz
	rz     = sum(r*z)
        p      = z + p * (rz/rz_old)
	
    
!         write(line,*) norm_res2
! 	k = mexPrintf(line//achar(13))

      end do
   

      call mexErrMsgTxt('CG didn''t converge in maxiter iterations')
      
 100  continue
      
      iter = i
      norm1 = sqrt(norm_res2)
      norm2 = maxres/(maxx-minx)
      
 110  continue
      
      deallocate(p)
      deallocate(b)
      deallocate(r)
      deallocate(r_old)
      deallocate(temp_vec)
      deallocate(z)
      
      return
      end

C-------------------------------------------------------------------------------
C  y = Ax
C-------------------------------------------------------------------------------
      subroutine matvec_prod(N, A, dia, nd, x, y)
      implicit none
      
      integer N, nd, dia(nd)

      double precision A(N,nd), x(N), y(N)

      integer i, j

C     main diagonal
      do i=1,N
        y(i) = A(i,1)*x(i)
      end do
      
      if (nd.gt.1) then
          do j=2,nd

            do i=1,N-dia(j)
              y(i)        = y(i)        + A(i,j)*x(i+dia(j))
              y(i+dia(j)) = y(i+dia(j)) + A(i,j)*x(i)
            end do

          end do
           
      end if
      
      return
      end

C-------------------------------------------------------------------------------
C  x: Mx=y
C-------------------------------------------------------------------------------
      subroutine precond(N, A, dia, nd, D, y, x)
      implicit none
      
      integer N, nd, dia(nd), dia2(nd+1)
      
      double precision A(N,nd), D(N), y(N)
      double precision x(N)
      
      double precision, dimension(:), allocatable  :: q
      double precision s
      
      integer i, j, k
      
      allocate(q(N))
            
      dia2(1:nd) = dia(1:nd)
      dia2(nd+1) = N

C--------------------
C  Diagonal scaling
C--------------------
C      do i=1,N
C        x(i) = y(i) / A(i,1)
C      end do
C-----------------------
C  Incomplete Choleski
C-----------------------
      i=1
      q(i) = y(i) / D(i)
        
      if (nd.gt.1) then
        do j=2,nd

            do i=dia2(j)+1,dia2(j+1)

                s = 0.d0
                do k=2,j
                    s = s - A(i-dia2(k),k)*q(i-dia2(k))
                end do
                q(i) = (y(i) + s ) / D(i)

            end do
                        
        end do
      
      end if 

      i=N
      x(i) = q(i)
      
      if (nd.gt.1) then
         
          do j=2,nd
             
              do i=N-dia2(j),N-dia2(j+1)+1,-1

                  s = 0.d0
                  do k=2,j
                     s = s - A(i,k)*x(i+dia2(k))                     
                  end do
                      
                  x(i) = q(i) + s / D(i)
                  
              end do
        
          end do
          
      end if

  
      deallocate(q)
      
      return
      end
C-------------------------------------------------------------------------------
C  
C-------------------------------------------------------------------------------
      subroutine make_cholinc(N, A, dia, nd, D)
      implicit none
      
      integer N, nd, dia(nd), dia2(nd+1)
      
      double precision A(N,nd)
      double precision D(N)
      
      integer i, j, k

      double precision s
      
      dia2(1:nd) = dia(1:nd)
      dia2(nd+1) = N

      i=1
      D(i) = A(i,1)
            
      do j=2,nd
         
          do i=dia2(j)+1,dia2(j+1)
          
              s = 0.d0

              do k=2,j
                  s = s - (A(i-dia2(k),k)**2.d0)/D(i-dia2(k));
              end do
                  
              D(i) = A(i,1) + s
              
          end do
      end do
      
C      write(6,"(E25.15)") (D(i), i=1,N)
  
      return
      end