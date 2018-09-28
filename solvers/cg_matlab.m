% cg_matlab
%-------------------------------------------------------------------------------
%  solve for x: Ax = y
%-------------------------------------------------------------------------------
% A should be sparse and spd

function [x,iter,norm1,norm2] = cg_matlab(A,y,tol,x_in,maxiter,D)

      % input:
      % A: input matrix, sparse Matlab format
      % y: rhs vector
      % tolerance
      % N
      % x_in: initial guess
      % max_iter: maximum number of iterations
      % D: (optional) preconditioning vector, created by make_cholinc
      
      % output
      % x: solution vector
      % iter: number of iterations
      % norm1 and norm2: error norms
      
      N     = size(A,1);
      N2    = size(A,2);
      
      if (N~=N2)
          error('matrix A should be square');
      end
      
      [B,d] = spdiags(A);
      nd    = (length(d)+1)/2;
      % take only lower diagonals
      B     = B(:,nd:-1:1);
      d     = d(nd:end);
      tol2  = tol*tol;

      if (nd==1 && d(1)==0) % only main diagonal
          x = y./B;
          iter = 1;
          norm1 =  norm(B.*x-y ) / norm(y);
          norm2 =  norm(B.*x-y,'inf') / norm(y,'inf');
          return
      end
      
%     fill x with x_in
      x  = x_in;
      
%     set up preconditioner if not provided
      if (nargin<9)
          D  = make_cholinc(A);
      end
      
      temp_vec =  matvec_prod(N, B, d, x);
      
      r  = y-temp_vec;
      z  = precond(N, B, d, D, r);     
      p  = z;
      norm_res2_0 = sum(y.*y);
      norm_res2   = norm_res2_0;

      if (norm_res2<=tol2)
        iter = 0;
        norm1 = 0.;
        norm2 = 0.;
        return
      end
      
      rz  = sum(r.*z);

      i = 1;
      while (norm_res2>tol2)
	
        temp_vec = matvec_prod(N, B, d, p);
        temp     = sum(p.*temp_vec);
        alpha    = rz/temp;
        
        x = x+alpha*p;
        r = r-alpha*temp_vec;
        
        norm_res2 = sum(r.^2);
        norm_res2 = norm_res2 / norm_res2_0;

        maxx = max(x);
        minx = min(x);
        maxres = max(abs(r));

        if (i>maxiter) 
            norm1 = sqrt(norm_res2);
            error(['CG didn''t converge in ' num2str(maxiter) ' iterations; residual= ' num2str(norm1)]);
        end
	
        z      = precond(N, B, d, D, r);
	    rz_old = rz;
        rz     = sum(r.*z);
	    p      = z+(rz/rz_old)*p;
        i      = i+1;
        
      end
              
      iter  = i-1;
      norm1 = sqrt(norm_res2);
      norm2 = maxres/(maxx-minx);

end


%-------------------------------------------------------------------------------
%  y = Ax
%-------------------------------------------------------------------------------
function y = matvec_prod(N, A, d, x)

      nd   = length(d);
      % A consists of lower diagonals of original matrix

     
      % main diagonal
      i = 1:N;
      y = A(i,1).*x(i);
      
      if (nd>1)
          for j=2:nd

            i         = 1:N-d(j);
            y(i)      = y(i)      + A(i,j).*x(i+d(j));
            y(i+d(j)) = y(i+d(j)) + A(i,j).*x(i); 

          end
           
      end
end

%-------------------------------------------------------------------------------
%  x: Mx=y
%-------------------------------------------------------------------------------
function x = precond(N, A, d, D, y)
      
     
      x    = zeros(N,1);
      q    = zeros(N,1);
      
      d    = [d; N];
      nd   = length(d);
      
%--------------------
%  Diagonal scaling
%--------------------
%       x   = y./A(:,1);

%-----------------------
%  Incomplete Choleski
%-----------------------
      i=1;
	  q(i) = y(i) / D(i);
      
      if (nd>2)
        
        for j=2:nd-1

            for i=d(j)+1:d(j+1)

                s=0;
                for k=2:j
                    s = s - A(i-d(k),k)*q(i-d(k));
                end
                q(i) = (y(i) + s ) / D(i);

            end
                        
        end
      
      end
 
      i=N;
      x(i) = q(i);
      
      if (nd>2)
         
          for j=2:nd-1
             
              for i=N-d(j):-1:N-d(j+1)+1
%                  
                  s=0;
                  for k=2:j
                     s = s - A(i,k)*x(i+d(k));                      
                  end
%                   k=2:j;
%                   s=- A(i,k)*x(i+d(k));
                      
                  x(i) = q(i) + s / D(i);
                  
              end
        
          end
          
      end

        
end
%-------------------------------------------------------------------------------
%  
%-------------------------------------------------------------------------------
