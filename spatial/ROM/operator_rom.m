function options = operator_rom(options)
% we construct the ROM operators in  'non-intrusive' way, i.e.
% we hardly require any knowledge of the exact FOM discretization. we
% simply make calls to diffusion and convection with the modes phi as input
% arguments

NV  = options.grid.Nu+options.grid.Nv;
B   = options.rom.B;

if (options.rom.weighted_norm == 0)
    Diag = options.grid.Om_inv;
elseif (options.rom.weighted_norm == 1)
    Diag = ones(NV,1);
end

% this is the projector for the momentum equations:
P = B'*spdiags(Diag,0,NV,NV);



%% diffusion
[yDiff,Diff] = operator_rom_diffusion(P,options);

options.rom.Diff  = Diff;
options.rom.yDiff = yDiff;

%% convection 
[conv_bc,conv_linear,conv_quad] = operator_rom_convection(P,options);

options.rom.Conv_quad   = conv_quad;
options.rom.Conv_linear = conv_linear;
options.rom.yConv       = conv_bc;

%% body force
% construct at t=t_start
[Fx, Fy] = force(options.time.t_start,options);
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
        [Fx, Fy] = force(0,options);
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
