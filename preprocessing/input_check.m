%% input checking
if (order4 == 1)
    if (~strcmp(visc,'laminar'))
        error('order 4 only implemented for laminar flow; need to add Su_vx and Sv_uy for 4th order method');
    end
    
    if (regularize~=0)
        error('order 4 only implemented for standard convection with regularize=0');
    end
end


%% initialize solution vectors

if (restart.load == 1)
    
    %% load restart file
    
    % load u,v,w,p,n,t,dt
    nzeros = filelen-length(num2str(restart.file));
    n_new  = num2str(10^nzeros);
    load ([restart.folder '/restart_' n_new(2:end) num2str(restart.file) '.mat']);
    
    % 'chop' off the part of convergence.txt and iterations_pressure.txt
    % that have iterations > n
    convergence_temp  = fileread([restart.folder '/convergence.txt']);
    pos               = regexp(convergence_temp,['\n' num2str(n+1)],'once');
    fconv_temp        = fopen([restart.folder '/convergence_temp.txt'],'w+');
    fwrite(fconv_temp,convergence_temp(1:pos));
    unix(['mv -f ' restart.folder '/convergence_temp.txt ' restart.folder '/convergence.txt']);
    
    iterations_temp   = fileread([restart.folder '/iterations_pressure.txt']);
    pos               = regexp(iterations_temp,['\n' num2str(n+1)],'once');
    fpres_temp        = fopen([restart.folder '/iterations_pressure_temp.txt'],'w+');
    fwrite(fpres_temp,iterations_temp(1:pos));
    unix(['mv -f ' restart.folder '/iterations_pressure_temp.txt ' restart.folder '/iterations_pressure.txt']);
    
else
    
    % loop index
    n     = 1;
    
    
    %% initial velocity field
%     uh    = u_start(:);
%     vh    = v_start(:);
    V     = V_start; %[uh; vh];
    V_old = V;
    
    
    switch visc
        case 'keps'
            kth = kt(:);
            eh  = e(:);
    end
    
    % for unsteady problems allocate k, umom and vmom, maxdiv and time
    % for steady problems the time for allocating during running is negligible,
    % since normally only a few iterations are required
    if (steady==0)
        t     = t_start;
        
        if (timestep.set==1)
            set_timestep;
        end
        % estimate number of time steps that will be taken
        nt     = ceil((t_end-t_start)/dt);
%         nsample = 
        % allocate variables, including initial condition
        maxres = zeros(nt+1,1);
        maxdiv = zeros(nt+1,1);
        k      = zeros(nt+1,1);
        umom   = zeros(nt+1,1);
        vmom   = zeros(nt+1,1);
        time   = zeros(nt+1,1);
        nonlinear_its = zeros(nt+1,1);
        
    end
    
    %% temperature
    switch options.case.boussinesq
        case 'temp'
            T = T_start(:);
        case 'none'
            T = 0;
        otherwise
            error('wrong setting for temperature');
    end   
    
    %% kinetic energy and momentum of initial velocity field
    % iteration 1 corresponds to t=0 (for unsteady simulations)
    [maxdiv(1), umom(1), vmom(1), k(1)] = check_conservation(V,T,t,options);
    
    
    if (maxdiv(1)>1e-12 && steady==0)
        fprintf(fcw,['warning: initial velocity field not (discretely) divergence free: ' num2str(maxdiv(1)) '\n']);
        fprintf(fcw,'additional projection to make initial velocity field divergence free\n');
        
        % make velocity field divergence free
        Om_inv = options.grid.Om_inv;
        G  = options.discretization.G;
        Nu = options.grid.Nu;
        Nv = options.grid.Nv;
        
        f  = options.discretization.M*V + options.discretization.yM;
        dp = pressure_poisson(f,t,options);
        V  = V - Om_inv.*(G*dp);
        %         uh = V(1:Nu); vh = V(Nu+1:end);
        % repeat conservation with updated velocity field
        [maxdiv(1), umom(1), vmom(1), k(1)] = check_conservation(V,T,t,options);
    end
    
    [symmetry_flag, symmetry_error] = check_symmetry(V,t,options);       
    
    
    %% initialize pressure
    if (options.rom.rom == 0)
        if (options.case.steady==1)
            % for steady state computations, the initial guess is the provided initial condition
            p  = p_start(:);
        else
            if (options.solversettings.p_initial == 1)
                % calculate initial pressure from a Poisson equation
                p = pressure_additional_solve(V,p_start(:),T,t,options);
            else
                % use provided initial condition (not recommended)
                p = p_start(:);
            end
        end
    end
    
    % ROM: uses the IC for the pressure; note that in solver_unsteady the pressure can be
    % computed from the ROM PPE, after the ROM basis has been set-up
    if (options.rom.rom==1)
        if (options.case.steady==1)
            error('ROM not implemented for steady flow');
        else
            p = p_start(:); 
        end
        
        if (options.BC.BC_unsteady==1 && options.rom.rom_bc < 2)
            error('unsteady BC typically requires rom_bc=2');
        end
        if (options.rom.rom_bc== 2 && options.BC.BC_unsteady==0)
            warning('unsteady BC set for ROM, but BC are steady');
        end
    end
    
    %%
    % for steady problems, with Newton linearization and full Jacobian,
    % first start with nPicard Picard iterations
    if (options.case.steady==1)
        options.solversettings.Newton_factor = 0;
    elseif (options.case.steady==0)
        % implicit RK time integration: set Newton_factor to 1
        if ( (method==21 || method==22) || (exist('method_startup','var') && method_startup==21)) 
            options.solversettings.Newton_factor = 1;
        end
    else
        error('wrong setting for steady parameter');
    end
    
    
    %% residual of momentum equations at start
%     if (strcmp(options.case.visc,'laminar'))
        [maxres(1), ~, ~] = F(V,V,p,T,t,options,0);
%     else
%         maxres(1) = 0;
%     end
    
    if (steady==0 && save_unsteady == 1)
        if (options.output.save_eco == 0)
%             % allocate space for variables
%             if (steady==0 && save_unsteady == 1)
%             ntemp = nt-sampling_start;    
%             uh_total = zeros(ntemp,options.grid.Nu);
%             vh_total = zeros(ntemp,options.grid.Nv);
%             p_total  = zeros(ntemp,options.grid.Np);
%             else
            uh_total = zeros(nt,options.grid.Nu);
            vh_total = zeros(nt,options.grid.Nv);
            p_total  = zeros(nt,options.grid.Np);
            T_total  = zeros(nt,options.grid.NT);
%             end                               
            % store initial solution
            V_start = V; % can have changed due to div-correction
            uh_total(1,:) = V_start(options.grid.indu);
            vh_total(1,:) = V_start(options.grid.indv);
            p_total(1,:)  = p_start(:);
            switch options.case.boussinesq   
                case 'temp'
                    T_total(1,:) = T_start(:);
            end
        
        elseif (options.output.save_eco == 1)
            mkdir([options.output.path_results '/snapshots']);
            save([options.output.path_results '/snapshots/Vp_n=1.mat'],'V','p');
        end

    end
                % Initialize statistics 
    switch options.case.boussinesq
        case 'temp'
            %             if (n>nsample && n==nsample_end)
            NV=options.grid.NV;
            Np=options.grid.Np;
            NT=options.grid.NT;
            Vmean=zeros(NV,1);
            pmean=zeros(Np,1);
            Tmean=zeros(NT,1);
            V_var=zeros(NV,1);
            p_var=zeros(Np,1);
            T_var=zeros(NT,1);         
%                     options.savestatisticsRBC.Vmean=Vmean;
%                     pmean=zeros(Np,1),options.savestatisticsRBC.Vmean=Vmean;
%                     Tmean=zeros(NT,1);
            end
    
end

if (rtp.show == 0 && rtp.movie==1)
    warning('real-time plotting off (rtp.show=0) but movie generation is on (rtp.movie=1). no movie will be generated.');
end

if (steady==0)
    if (save_unsteady==1 && save_file==0)
        warning('unsteady data is stored in workspace (save_unsteady=1) but not saved to a file (save_file=0)');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
