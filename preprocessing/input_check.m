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
    uh    = u_start(:);
    vh    = v_start(:);
    V     = [uh; vh];
    V_old = V;
    
    
    if (strcmp(visc,'turbulent'))
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
        % estimate number of time steps
        nt     = ceil((t_end-t)/dt);
        
        maxres = zeros(nt,1);
        maxdiv = zeros(nt,1);
        k      = zeros(nt,1);
        umom   = zeros(nt,1);
        vmom   = zeros(nt,1);
        time   = zeros(nt,1);
        
    end
    
    
    %% kinetic energy and momentum of initial velocity field
    % iteration 1 corresponds to t=0 (for unsteady simulations)
    [maxdiv(1), umom(1), vmom(1), k(1)] = check_conservation(V,t,options);
    
    
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
        [maxdiv(1), umom(1), vmom(1), k(1)] = check_conservation(V,t,options);
    else
        [symmetry_flag, symmetry_error] = check_symmetry(uh,vh,t,options);
        
    end
    
    
    %% body forces
    
    % set body force if not set
    %     if (exist('Fx','var')~=1)
    %         Fx    = zeros(Nu,1);
    %     end
    %     if (exist('Fy','var')~=1)
    %         Fy    = zeros(Nv,1);
    %     end
    
    %     if (strcmp(options.case.visc,'turbulent'))
    %         % initialize body forces k and epsilon equation
    %         Fk  = zeros(Np,1);
    %         Fe  = zeros(Np,1);
    %     end
    
    
    %% initialize pressure
    if (steady==1 || rom==1)
        % for steady state computations, the initial guess is the provided initial condition
        p  = p_start(:);
    else
        if (p_initial == 1)
            % calculate initial pressure from a Poisson equation
            p = pressure_additional_solve(V,p_start(:),t,options);
        else
            % use provided initial condition (not recommended)
            p = p_start(:);
        end
    end
    
    %%
    % for steady problems, with Newton linearization and full Jacobian,
    % first start with nPicard Picard iterations
    if (steady==1)
        options.solversettings.Newton_factor = 0;
    elseif (steady==0)
        if (method==21 || method_startup==21) % implicit RK time integration
            options.solversettings.Newton_factor = 1;
        end
    end
    
    %% residual of momentum equations at start
    if (strcmp(options.case.visc,'laminar'))
        [maxres(1), ~, ~] = F(V,p,t,options,0);
    else
        maxres(1) = 0;
    end
    
    if (steady==0 && save_unsteady == 1)
        % allocate space for variables
        uh_total = zeros(nt,options.grid.Nu);
        vh_total = zeros(nt,options.grid.Nv);
        p_total  = zeros(nt,options.grid.Np);
        % store initial solution
        uh_total(1,:) = uh;
        vh_total(1,:) = vh;
        p_total(1,:)  = p;
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