%  Main solver file for unsteady calculations

%% load restart file if necessary
if (restart.load == 0 && options.output.save_results==1)
    fprintf(fconv,'n            dt               t                res              maxdiv           umom             vmom             k\n');
    fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n',...
        n,dt,t,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
end

%% plot initial solution

% compute V_bc snapshot if BC unsteady and snapshots are wanted
Om_inv = options.grid.Om_inv;
if (options.BC.BC_unsteady == 1 && options.rom.rom_bc == 2 && options.rom.pro_rom == 1)
%     options = set_bc_vectors(t,options);
%     f       = options.discretization.yM;
%     dp      = pressure_poisson(f,t,options);
%     Vbc(:,n)= - Om_inv.*(options.discretization.G*dp);
Vbc(:,n) = get_unsteadyVbc(t,options);
end


if (rtp.show==1)
    run(rtp.file);
    % for movies, capture this frame
    if (rtp.movie==1)
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
end

%% start-up of multistep methods

% for methods that need convection from previous time step
if (method==2)
    if (options.BC.BC_unsteady == 1)
        options = set_bc_vectors(t,options);
    end
    [convu_old,convv_old] = convection(V,V,t,options,0);
    conv_old = [convu_old;convv_old];
end

% for methods that need u^(n-1)
if (method == 5)
    V_old  = V;
    p_old  = p;
end

% set current velocity and pressure
Vn = V;
pn = p;
tn = t;

% for methods that need extrapolation of convective terms
if (method == 62 || method == 92 || method==142 || method==172 || method==182 || method==192)
    V_ep      = zeros(Nu+Nv,method_startup_no);
    V_ep(:,1) = V;
end


dtn    = dt;

% eps    = 1e-12;

method_temp = method;

% addpath 'C:\Users\20201213\Documents\Uni\Master thesis\clean_code\INS2D\debug_stuff'
% addpath 'debug_stuff'
% jacobian_test


%% start time stepping
time_start = toc

k_sum2 = zeros(nt,1);

while(n<=nt)
    
    % time step counter
    n = n+1;
    
    % for methods that need a velocity field at n-1 the first time step
    % (e.g. AB-CN, oneleg beta) use ERK or IRK 

    if ((method_temp==2 || method_temp==5) && n<=method_startup_no)
        fprintf(fcw,['starting up with method ' num2str(method_startup) '\n']);
        method      = method_startup;        
    else
        method      = method_temp;
    end
    
    % time reversal (used in inviscid shear-layer testcase)
%     if (abs(t-t_end/2)<1e-12)
%         dt = -dt;
%     end
%         
    
    % perform a single time step with the time integration method
    if (method==2)
        [V,p,conv] = time_AB_CN(Vn,pn,conv_old,tn,dt,options);
        conv_old = conv;
    elseif (method==5)
        [V,p] = time_oneleg(Vn,pn,V_old,p_old,tn,dt,options);

    elseif (method==20)
        [V,p] = time_ERK(Vn,pn,tn,dt,options);
    elseif (method==21)
        [V,p,nonlinear_its(n),k_sum2(n)] = time_IRK(Vn,pn,tn,dt,options);
    else
        error('time integration method unknown');
    end
    
    
    % the velocities and pressure that are just computed are at
    % the new time level t+dt:
    t = tn + dt;
    time(n) = t;
    
    % update solution
    V_old = Vn;
    p_old = pn;                
    Vn = V;
    pn = p;
    tn = t;      
    
    % compute V_bc snapshot if BC unsteady and snapshots are wanted
    if (options.BC.BC_unsteady == 1 && options.rom.pro_rom == 1)
        options = set_bc_vectors(t,options);
        f       = options.discretization.yM;
        dp      = pressure_poisson(f,t,options);
        Vbc(:,n)= - Om_inv.*(options.discretization.G*dp);
    end
    
    % check residuals, conservation, set timestep, write output files
    process_iteration;
    
    
end
disp('finished time-stepping...');
time_loop = toc-time_start