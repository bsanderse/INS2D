function options = operator_postprocessing(options)
% construct postprocessing operators such as vorticity

% boundary conditions
BC = options.BC;

Nx = options.grid.Nx;
Ny = options.grid.Ny;

hx = options.grid.hx;
hy = options.grid.hy;
gx = options.grid.gx;
gy = options.grid.gy;

gxi = options.grid.gxi;
gyi = options.grid.gyi;
gxd = options.grid.gxd;
gyd = options.grid.gyd;


order4 = options.discretization.order4;

if (order4==1)
    alfa = options.discretization.alfa;
    beta = options.discretization.beta;
    gxi3 = options.grid.gxi3;
    gyi3 = options.grid.gyi3;
    Omvort3 = options.grid.Omvort3;
end

%% vorticity

% operators act on internal points only
%
% du/dy, like Su_uy
diag1  = 1./gy(2:end-1);
W1D    = spdiags([-diag1 diag1],[0 1],Ny-1,Ny);
% extend to 2D
Wu_uy  = kron(W1D,speye(Nx-1));

% dv/dx, like Sv_vx
diag1  = 1./gx(2:end-1);
W1D    = spdiags([-diag1 diag1],[0 1],Nx-1,Nx);
% extend to 2D
Wv_vx  = kron(speye(Ny-1),W1D);


%% for periodic BC, covering entire mesh
if (strcmp(BC.u.left,'per') && strcmp(BC.v.low,'per'))

    % du/dy, like Su_uy
    diag1  = 1./gyd;
    W1D    = spdiags([diag1 -diag1 diag1 -diag1],[-Ny -1 0 Ny-1],Ny+1,Ny);
    % extend to 2D
    diag2  = ones(Nx,1);
    Wu_uy  = kron(W1D,spdiags([diag2 diag2],[-Nx 0],Nx+1,Nx));

    % dv/dx, like Sv_vx
    diag1  = 1./gxd;
    W1D    = spdiags([diag1 -diag1 diag1 -diag1],[-Nx -1 0 Nx-1],Nx+1,Nx);
    % extend to 2D
    diag2  = ones(Ny,1);
    Wv_vx  = kron(spdiags([diag2 diag2],[-Ny 0],Ny+1,Ny),W1D);

end


% %% for periodic BC, such that sum(Omega_w)=L_x*L_y (for studying conservation of omega)
% if (strcmp(BC.u.left,'per') && strcmp(BC.v.low,'per'))
%     Omvort = options.grid.Omvort;
%     % du/dy, like Su_uy
%     diag1  = 1./gyi;
%     W1D    = spdiags([-diag1 diag1],[0 1],Ny,Ny+1);
%     W1D_BC = BC_general_stag(Ny+1,Ny,1,BC.u.low,BC.u.up,hy(1),hy(end));
%     % extend to 2D
%     % point value
%     Wu_uy_p= kron(W1D*W1D_BC.B1D,speye(Nx));
%     Wu_uy  = spdiags(Omvort,0,Nx*Ny,Nx*Ny)*Wu_uy_p;
%     
%     % dv/dx, like Sv_vx
%     diag1  = 1./gxi;
%     W1D    = spdiags([-diag1 diag1],[0 1],Nx,Nx+1);
%     W1D_BC = BC_general_stag(Nx+1,Nx,1,BC.v.left,BC.v.right,hx(1),hx(end));
%     % extend to 2D
%     % point value
%     Wv_vx_p= kron(speye(Ny),W1D*W1D_BC.B1D);
%     Wv_vx  = spdiags(Omvort,0,Nx*Ny,Nx*Ny)*Wv_vx_p;
%     
%     if (order4==1)
%         %     diag1  = 1./gyd3;
%         %     W1D    = spdiags([diag1 -diag1 diag1 -diag1],[-Ny+1 -2 1 Ny-2],Ny,Ny);
%         %     % extend to 2D
%         %     Wu_uy3 = kron(W1D,speye(Nx));
%         diag1  = 1./gyi3;
%         W1D    = spdiags([-diag1 diag1],[0 3],Ny,Ny+3);
%         W1D_BC = BC_vort3(Ny+3,Ny,3,BC.u.low,BC.u.up,hy(1),hy(end));
%         
%         % extend to 2D
%         Wu_uy3 = kron(W1D*W1D_BC.B1D,speye(Nx));
%         
%         % dv/dx, like Sv_vx
%         diag1  = 1./gxi3;
%         W1D    = spdiags([-diag1 diag1],[0 3],Nx,Nx+3);
%         W1D_BC = BC_vort3(Nx+3,Nx,3,BC.v.left,BC.v.right,hx(1),hx(end));
%         
%         % extend to 2D
%         Wv_vx3 = kron(speye(Ny),W1D*W1D_BC.B1D);
%         
%         % this is the conserved quantity:
%         Wu_uy1 = Wu_uy;
%         Wv_vx1 = Wv_vx;
%         Wu_uy  = (alfa*Wu_uy1  - spdiags(Omvort3,0,Nx*Ny,Nx*Ny)*Wu_uy3);
%         Wv_vx  = (alfa*Wv_vx1  - spdiags(Omvort3,0,Nx*Ny,Nx*Ny)*Wv_vx3);
%         
%         % this is a fourth order approximation to a point value:
%         Wu_uy_p  = beta*Wu_uy_p  + (1-beta)*Wu_uy3;
%         Wv_vx_p  = beta*Wv_vx_p  + (1-beta)*Wv_vx3;
%         
%         % this is a fourth order approximation to a volume integrated quantity
%         diag1   = ones(Nx,1);
%         W1Dx    = (1/24)*spdiags([diag1 10*diag1 diag1],[0 1 2],Nx,Nx+2);
%         W1D_BCx = BC_general(Nx+2,Nx,2,BC.u.left,BC.u.right,hx(1),hx(end));
%         diag1   = ones(Ny,1);
%         W1Dy    = (1/24)*spdiags([diag1 10*diag1 diag1],[0 1 2],Ny,Ny+2);
%         W1D_BCy = BC_general(Ny+2,Ny,2,BC.v.low,BC.v.up,hy(1),hy(end));
%         
%         
%         Wint    = kron(speye(Ny),W1Dx*W1D_BCx.B1D) + kron(W1Dy*W1D_BCy.B1D,speye(Nx));
%     end
%     
% end

options.discretization.Wv_vx = Wv_vx;
options.discretization.Wu_uy = Wu_uy;


end
