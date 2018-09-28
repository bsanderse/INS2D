function D = make_cholinc(A)
      
      N     = size(A,1);
      N2    = size(A,2);
      
      if (N~=N2)
          error('matrix A should be square');
      end

      [A,d] = spdiags(A);
      nd    = (length(d)+1)/2;
      % take only lower diagonals
      A     = A(:,nd:-1:1);
      d     = d(nd:end);

     
      d    = [d; N];
      nd   = length(d);
      D    = zeros(N,1);
      
      i=1;
      D(i) = A(i,1);
      
      
      for j=2:nd-1
         
          for i=d(j)+1:d(j+1)
          
              s = 0;

              for k=2:j
                  s = s - A(i-d(k),k)^2/D(i-d(k));
              end
                  
              D(i) = A(i,1) + s; 
              
          end
      end
%       for i=2:Nx
%         D(i) = A(i,1) - A(i-1,2)^2/D(i-1);
%       end
%       
%       for i=Nx+1:NxNy
%         D(i) = A(i,1) - A(i-1,2)^2/D(i-1) - A(i-Nx,3)^2/D(i-Nx);
%       end
%       
%       for i=NxNy+1:Nt
%         D(i) = A(i,1) - A(i-1,2)^2/D(i-1) - A(i-Nx,3)^2/D(i-Nx) ...
%                                           - A(i-NxNy,4)^2/D(i-NxNy);
%       end
     
end