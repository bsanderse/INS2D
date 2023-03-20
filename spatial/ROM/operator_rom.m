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

if (options.rom.precompute_convection == 1 || options.rom.precompute_diffusion == 1 || ...
    options.rom.precompute_force == 1)
    % this is the projector for the momentum equations:
    P = B'*spdiags(Diag,0,NV,NV);

end
% added for temperature 

switch options.case.boussinesq   
    case 'temp'
        BT  = options.rom.BT;
        NT  = options.grid.NT;
        MT  = options.rom.MT;
        if (options.rom.weighted_norm_T == 0)
            DiagT = options.grid.Omp_inv;
        elseif (options.rom.weighted_norm_T == 1)
            DiagT = ones(NT,1);
        end
        if (options.rom.precompute_convectionT == 1 || options.rom.precompute_diffusionT == 1)
            % this is the projector for the temperature equation:
            PT = BT'*spdiags(DiagT,0,NT,NT);
        end
end

%% diffusion
% if (options.rom.precompute_diffusion == 1 || options.rom.precompute_diffusionT == 1)
if (options.rom.precompute_diffusion == 1)
    [yDiff,Diff] = operator_rom_diffusion(P,options);
%     [yDiff,Diff,yDiffT,DiffT] = operator_rom_diffusion(P,PT,options);
   
    options.rom.Diff  = Diff;
    options.rom.yDiff = yDiff;
%     options.rom.DiffT = DiffT;
end

%% convection 
if (options.rom.precompute_convection == 1)
    [conv_bc,conv_linear,conv_quad] = operator_rom_convection(P,options);

    options.rom.Conv_quad   = conv_quad;
    options.rom.Conv_linear = conv_linear;
    options.rom.yConv       = conv_bc;
end


%% body force
% always precomputed if forcing is steady
if (options.force.isforce == 1)
    if (options.rom.precompute_force == 1)
        % construct at t=t_start with dummy velocity field
        [Fx, Fy] = force(zeros(NV,1),options.time.t_start,options,0);
        F        = P*[Fx;Fy];
        options.rom.F = F;
    end
else
    M  = options.rom.M;
    options.rom.F = zeros(M,1);
end

%% Include buoyancy force in the momentum equation
switch options.case.boussinesq
    
    case 'temp'
        % get T at v-locations
        % note that AT_v includes the volumes Omega_v
        Nu = options.grid.Nu;
        if (options.rom.precompute_buoyancy_force == 1)
            F_buoyancy_v_precompute = options.discretization.AT_v*BT;
            F_buoyancy_precompute = [zeros(Nu,MT); F_buoyancy_v_precompute];
            F_buoyancy_ROM_precompute = B'*(Diag.*F_buoyancy_precompute);
            options.rom.F_buoyancy_ROM_precompute = F_buoyancy_ROM_precompute;
        end
end
    
% FV = [Fu;Fv];

%% diffusion for temperature equation
if (options.rom.precompute_diffusionT == 1)
    [yDiffT,DiffT] = operator_rom_diffusionT(PT,options);
%     [DiffT] = operator_rom_diffusionT(PT,options);

    options.rom.DiffT  = DiffT;
    options.rom.yDiffT = yDiffT;
end

%% pressure
% the pressure gradient term in the momentum equation disappears in the ROM
% however, we still want to get the pressure, which is obtained by solving
% a Poisson equation on the ROM level
% the right hand side of this pressure equation consists of the ROM
% momentum equation projected onto the pressure basis
if (options.rom.div_free == 1)

    if (options.rom.pressure_recovery == 1)

        % generate Poisson matrix on ROM level
        Bp = options.rom.Bp;
        A_ROM = Bp'*options.discretization.A*Bp;
        % get LU decomposition
    %     [L,U] = lu(A_ROM);
    %     options.rom.L = L;
    %     options.rom.U = U;
        options.rom.A_decomp = decomposition(A_ROM);

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
    
elseif (options.rom.div_free == 0)
    
    % here we always precompute:
    options.rom.Mdiv = options.rom.Bp'*options.discretization.M*options.rom.B;
    options.rom.G = -options.rom.Mdiv';
    % note that the following Poisson matrix is different from the one in
    % pressure_additional solve (which only uses Bp for projection)
    A_ROM = options.rom.Mdiv*options.rom.G;
    options.rom.A_decomp = decomposition(A_ROM);
    
    % options.rom.yM contains the BC of the divergence equation for all
    % time steps
    options.rom.yMt = options.rom.Bp'*options.rom.yM;
    
    
end
    
