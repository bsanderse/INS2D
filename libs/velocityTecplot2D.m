function velocityTecplot2D(filename,x,y,u,v,p,T)
% function velocityTecplot2D(x,y,u,v,p,T,filename)
% input: u,v,p and T at x,y positions
% input: filename without extension

% output: tecplot file (.plt) with filename
   

    % Get the matrix dimensions.
    [Nx Ny] = size(u);
    Nz = 1;
    
    if (nargin<5)
        error('Missing input arguments');
    end
    
    if (nargin==5)
        % no temperature, no pressure
        T = zeros(Nx,Ny);
        p = zeros(Nx,Ny);
    elseif (nargin==6)
        % only pressure
        T = zeros(Nx,Ny);
    end
    
    
    % Open the file.
    fid = fopen([filename '.plt'], 'w');
    if fid == -1
        error('Cannot open file for writing.');
    end

    % New line.
    nl = sprintf('\n');

    % Write the file header.
    fwrite(fid, ['TITLE = "' num2str(filename) '"' nl ...
                 'VARIABLES  = "X", "Y", "U", "V", "P", "T"' nl...
                 'ZONE ZONETYPE=ORDERED, DATAPACKING=POINT, ' ...
                 'I=' num2str(Nx) ', J=' num2str(Ny) ', K=' num2str(Nz) nl]);
   
    xnew   = kron(ones(Ny,1),x);
    ynew   = kron(y,ones(Nx,1));
    varnew = [xnew(:) ynew(:) u(:) v(:) p(:) T(:)];    
    % datapacking=point:
    fprintf(fid, '%8e\t %8e\t %8e\t %8e\t %8e\t %8e\n', varnew');
  
%         for j=1:Ny
%             for i=1:Nx
%         
% %                 fwrite(fid,[num2str(x(i)) ' ' num2str(y(j)) ' ' num2str(z(k)) ' '...
% %                             num2str(u(i,j,k)) ' ' num2str(v(i,j,k)) ' ' num2str(w(i,j,k)) nl]);
%                 fprintf(fid, '%8e\t %8e\t %8e\t %8e\t %8e\t %8e\t %8e\t %8e\n', ...
%                              [x(i) y(j) 0 u(i,j) v(i,j) 0 p(i,j) T(i,j)]);
%                         
%             end
%         end

    fclose(fid);
end