function options = operator_rom(options)

NV  = options.grid.Nu+options.grid.Nv;
B   = options.rom.B;

if (options.rom.weighted_norm == 0)
    Diag = options.grid.Om_inv;
elseif (options.rom.weighted_norm == 1)
    Diag = ones(NV,1);
end

% this is the projector for the momentum equations:
P = B'*spdiags(Diag,0,NV,NV);

options.rom.P = P;

%% diffusion
% if options.rom.rom_bc == 2
    if options.rom.bc_recon == 0 
            [yDiff,Diff] = operator_rom_diffusion(P,options);
    
            options.rom.Diff  = Diff;
            options.rom.yDiff = yDiff;
    elseif options.rom.bc_recon == 1
        [yDiff,Diff,DiffBC] = operator_rom_diffusion_unsteadyBC(P,options);

        options.rom.Diff   = Diff;
        options.rom.DiffBC = DiffBC;
        options.rom.yDiff  = yDiff;
    else
            [yDiff2,Diff2,DiffBC2] = operator_rom_diffusion_unsteadyBC2(P,options);
  
            options.rom.Diff2   = Diff2;
            options.rom.DiffBC2 = DiffBC2;
            options.rom.yDiff2  = yDiff2;
    end
% else
%     [yDiff,Diff] = operator_rom_diffusion(P,options);
%     
%     options.rom.Diff  = Diff;
%     options.rom.yDiff = yDiff;
% end
%% convection 
% if options.rom.rom_bc == 2
    if options.rom.bc_recon == 0 
        [conv_bc,conv_linear,conv_quad] = operator_rom_convection(P,options);

        options.rom.Conv_quad   = conv_quad;
        options.rom.Conv_linear = conv_linear;
        options.rom.yConv       = conv_bc;
    elseif options.rom.bc_recon == 1
        [C_hom,C_hom_inhom,C_hom_bc,C_inhom,C_inhom_bc,C_bc] = ...
            operator_rom_convection_unsteadyBC(P,options);

        options.rom.C_hom       = C_hom;
        options.rom.C_hom_inhom = C_hom_inhom;
        options.rom.C_hom_bc    = C_hom_bc;
        options.rom.C_inhom     = C_inhom;
        options.rom.C_inhom_bc  = C_inhom_bc;
        options.rom.C_bc        = C_bc;
    else
        [C_hom2,C_hom_inhom2,C_hom_bc2,C_inhom2,C_inhom_bc2,C_bc2] = ...
            operator_rom_convection_unsteadyBC2(P,options);
        
        options.rom.C_hom2       = C_hom2;
        options.rom.C_hom_inhom2 = C_hom_inhom2;
        options.rom.C_hom_bc2    = C_hom_bc2;
        options.rom.C_inhom2     = C_inhom2;
        options.rom.C_inhom_bc2  = C_inhom_bc2;
        options.rom.C_bc2        = C_bc2;
    end
% else
%     [conv_bc,conv_linear,conv_quad] = operator_rom_convection(P,options);
%     
%     options.rom.Conv_quad   = conv_quad;
%     options.rom.Conv_linear = conv_linear;
%     options.rom.yConv       = conv_bc;
% end


%% body force
% construct at t=t_start with dummy velocity field
[Fx, Fy] = force(zeros(NV,1),options.time.t_start,options,0);
F        = P*[Fx;Fy];
options.rom.F = F;

%% pressure
% the pressure gradient term in the momentum equation disappears in the ROM
% however, we still want to get the pressure, which is obtained by solving
% a Poisson equation on the ROM level
% the right hand side of this pressure equation consists of the ROM
% momentum equation projected onto the pressure basis
if (options.rom.pressure_recovery == 1)

    % generate Poisson matrix on ROM level
    Bp = options.rom.Bp;
    A_ROM = Bp'*options.discretization.A*Bp;
    % get LU decomposition
    [L,U] = lu(A_ROM);
    options.rom.L = L;
    options.rom.U = U;
        
    
    if (options.rom.pressure_precompute == 1)
        % operators for right-hand side pressure equation

        % this is the projector for the Poisson equation:
        P_PPE = Bp'*options.discretization.M * spdiags(options.grid.Om_inv,0,NV,NV);    

        [conv_bc,conv_linear,conv_quad] = operator_rom_convection(P_PPE,options);
        [yDiff,Diff] = operator_rom_diffusion(P_PPE,options);    
        [Fx, Fy] = force(zeros(NV,1),options.time.t_start,options,0);
        F = P_PPE*[Fx;Fy];

        % note: for sign convention see F.m or F_ROM.m
        
        % constant terms in rhs 
        % we distinguish between force and BC, in order to allow time
        % varying forces
        options.rom.ppe_force  =  F;
        options.rom.ppe_bc     = -conv_bc + yDiff;
        % terms to be multiplied with R
        options.rom.ppe_linear = -conv_linear + Diff;
        % terms to be multiplied with kron(R,R)
        options.rom.ppe_quad   = -conv_quad;   
        
        % this is useful to do projection of time-varying quantities, e.g.
        % the force
        options.rom.P_PPE = P_PPE;
        
        % this is useful for evaluating int ( p u*n ) dS (pressure work):
        options.rom.yM = Bp'*options.discretization.yM;
    end
    
    if (options.rom.pressure_mean == 1)
        options.rom.ppe_mean = Bp'*options.discretization.A*options.rom.p_mean;
    end

end

%% open boundary condition contribution
% check whether open boundary conditions occur
BC = options.BC;
obc = max(strcmp({BC.u.left,BC.u.right,BC.u.low,BC.u.up},'mvp-obc')); 
if obc
%     if options.rom.rom_bc == 2
%         if options.rom.bc_recon == 3

%             [obc_hom,obc_inhom,obc_hom_inhom2,obc_hom_inhom,obc_inhom_hom] ...
            [obc_hom,obc_inhom,obc_hom_inhom2] ...
            = operator_rom_obc(P,options);
        
            options.rom.obc_hom = obc_hom;
            options.rom.obc_inhom = obc_inhom;
            options.rom.obc_hom_inhom2 = obc_hom_inhom2;
%             options.rom.obc_hom_inhom = obc_hom_inhom;
%             options.rom.obc_inhom_hom = obc_inhom_hom;
            
            switch BC.gO_type
                case 0
                    % nothin to do
                case 1
                    [gO_hom,gO_inhom] = operator_rom_gO(P,options);
                    
                    options.rom.gO_hom   = gO_hom;
                    options.rom.gO_inhom = gO_inhom;
                case 2
                    error('Sorry, obc offline decomposition for more complex gO not implemented')
            end
%         else
%             error('Sorry, precomputation not implemented for bc_recon =/= 3')
%         end
%     else
%         error('Sorry, precomputation not implemented for rom_bc =/= 2')
%     end
end

% compute matrix F_M which maps ybc to yM
% for 
