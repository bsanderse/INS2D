%% Divergence operator M

% note that the divergence matrix M is not square
mat_hx = spdiags(hx,0,Nx,Nx);
mat_hy = spdiags(hy,0,Ny,Ny);

% for fourth order: mat_hx3 is defined in operator_interpolation

%% Mx
% building blocks consisting of diagonal matrices where the diagonal is
% equal to constant per block (hy(block)) and changing for next block to
% hy(block+1)
diag1       = ones(Nux_t,1);
M1D         = spdiags([-diag1 diag1],[0 1],Nux_t-1,Nux_t);
                                  
% we only need derivative at inner pressure points, so we map the resulting
% boundary matrix (restrict)
diagpos = 0;
if (strcmp(BC.u.right,'pres') && strcmp(BC.u.left,'pres') )
    diagpos = 1;
end
if (~strcmp(BC.u.right,'pres') && strcmp(BC.u.left,'pres') )
    diagpos = 1;
end
if (strcmp(BC.u.right,'pres') && ~strcmp(BC.u.left,'pres') )
    diagpos = 0;
end
if (strcmp(BC.u.right,'per') && strcmp(BC.u.left,'per') ) % like pressure left
    diagpos = 1;
end

BMx    = spdiags(ones(Npx,1),diagpos,Npx,Nux_t-1);
M1D    = BMx*M1D;

% extension to 2D to be used in post-processing files
Bup    = kron(speye(Nuy_in),BMx);

% boundary conditions
Mx_BC  = BC_general(Nux_t,Nux_in,Nux_b, ...
                                      BC.u.left,BC.u.right,hx(1),hx(end)); 
Mx_BC.Bbc = kron(mat_hy,M1D*Mx_BC.Btemp);
                                  
% extend to 2D
Mx     = kron(mat_hy,M1D*Mx_BC.B1D);

if (order4==1)
    diag1       = ones(Nux_t+1,1);
    M1D3        = spdiags([-diag1 diag1],[0 3],Nux_t-1,Nux_t+2);
    M1D3        = BMx*M1D3;
    Mx_BC3      = BC_div2(Nux_t+2,Nux_in,Nux_t+2-Nux_in, ...
                                      BC.u.left,BC.u.right,hx(1),hx(end));
    Mx3         = kron(mat_hy3,M1D3*Mx_BC3.B1D);
    Mx_BC3.Bbc  = kron(mat_hy3,M1D3*Mx_BC3.Btemp);

end

%% My
% same as Mx but reversing indices and kron arguments
diag1  = ones(Nvy_t,1);
M1D    = spdiags([-diag1 diag1],[0 1],Nvy_t-1,Nvy_t);

% we only need derivative at inner pressure points, so we map the resulting
% boundary matrix (restriction)
diagpos = 0;
if (strcmp(BC.v.up,'pres') && strcmp(BC.v.low,'pres') )
    diagpos = 1;
end
if (~strcmp(BC.v.up,'pres') && strcmp(BC.v.low,'pres') )
    diagpos = 1;
end
if (strcmp(BC.v.up,'pres') && ~strcmp(BC.v.low,'pres') )
    diagpos = 0;
end
if (strcmp(BC.v.up,'per') && strcmp(BC.v.low,'per') ) % like pressure low
    diagpos = 1;
end

BMy    = spdiags(ones(Npy,1),diagpos,Npy,Nvy_t-1);
M1D    = BMy*M1D;
% extension to 2D to be used in post-processing files
Bvp    = kron(BMy,speye(Nvx_in));

% boundary conditions
My_BC  = BC_general(Nvy_t,Nvy_in,Nvy_b, ...
                                      BC.v.low,BC.v.up,hy(1),hy(end));
My_BC.Bbc = kron(M1D*My_BC.Btemp,mat_hx);

% extend to 2D
My     = kron(M1D*My_BC.B1D,mat_hx);

if (order4==1)
    diag1       = ones(Nvy_t+1,1);
    M1D3        = spdiags([-diag1 diag1],[0 3],Nvy_t-1,Nvy_t+2);
    M1D3        = BMy*M1D3;
    My_BC3      = BC_div2(Nvy_t+2,Nvy_in,Nvy_t+2-Nvy_in, ...
                                      BC.v.low,BC.v.up,hy(1),hy(end));
    My3         = kron(M1D3*My_BC3.B1D,mat_hx3);
    My_BC3.Bbc  = kron(M1D3*My_BC3.Btemp,mat_hx3);

end

    
%% resulting divergence matrix
if (order4==1)
    Mx     = alfa*Mx-Mx3;
    My     = alfa*My-My3;
end
M      = [Mx My];


%% Gradient operator G

% like in the continuous case, grad = -div*
Gx     = -Mx';
Gy     = -My';

G      = [Gx; Gy];

% if (ibm==1)
%     solver_unsteady_ibm;
% end


%% Pressure matrix for pressure correction method;
% also used to make initial data divergence free or compute additional poisson solve

if (steady == 0 && ~strcmp(visc,'turbulent'))
    
    %   Note that the matrix for the pressure is constant in time.
    %   Only the right hand side vector changes, so the pressure matrix 
    %   can be set up outside the time-stepping-loop.

    %   Laplace = div grad
    A     = M*spdiags(Om_inv,0,Nu+Nv,Nu+Nv)*G;

   %   LU decomposition
    if (poisson==3)
        if (exist(['cg.' mexext],'file')==3)
            [B, d]   = spdiags(A);
            ndia     = (length(d)+1)/2;
            dia      = d(ndia:end);
            B        = B(:,ndia:-1:1);               
        else
            fprintf(fcw,'No correct CG mex file available, switching to Matlab implementation \n');
            poisson = 4;
        end
    end
    if (poisson==4)
        % preconditioner
        A_pc = make_cholinc(A);
    end
    if (poisson==2)
        fprintf(fcw,'incomplete LU decomposition of pressure matrix...\n');
        setup.type = 'nofill';
        [L, U] = ilu(-A,setup);
    end
    if (poisson==1)
        fprintf(fcw,'LU decomposition of pressure matrix...\n');
        [L, U] = lu(A);
    end
    
    if (poisson==5)
        fprintf(fcw,'Petsc for pressure matrix solution... waiting for matrix\n');
        % open socket only once
%         system('petscmpiexec -n 2 ./solvers/petsc_poisson_par -viewer_socket_port 5600 -pc_type hypre -pc_hypre_type boomeramg &');
        PS = PetscOpenSocket(5600);
        PetscBinaryWrite(PS,-A);
        fprintf(fcw,'...matrix written to Petsc!\n');
    end

    % check if all the row sums of the pressure matrix are zero, which
    % should be the case if there are no pressure boundary conditions
    if (~strcmp(BC.v.low,'pres') && ~strcmp(BC.v.up,'pres') && ...
        ~strcmp(BC.u.right,'pres') && ~strcmp(BC.u.left,'pres'))    
        if (max(abs(A*ones(Np,1)))>1e-10)
            fprintf(fcw,'warning: pressure matrix: not all rowsums are zero!\n');
        end
    end

end