%  Main solver file for unsteady calculations

%% load restart file if necessary
if (restart.load == 0 && options.output.save_results==1)
    fprintf(fconv,'n            dt               t                res              maxdiv           umom             vmom             k\n');
    fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n',...
        n,dt,t,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
end

%% plot initial solution
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
V_old  = V;
p_old  = p;

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


%% start time stepping
% while (abs(t)<=(t_end-dt+eps))
% rev = 0;
while(n<=nt)
    
    % time reversal
    %       if (n==nt/2+1)
    %           fprintf(fcw,['reversing time at t=' num2str(t) '\n']);
    %           dt=-dt;
    %       end
    % time reversal for linearized methods
    %       if (n==nt/2+2 && rev==0)
    %           n = n-1;
    %           dt=-dt;
    %           t = t+dt;
    %           uh = uhn;
    %           vh = vhn;
    %           V  = [uh;vh];
    %           p  = pn;
    %           fprintf(fcw,['reversing time at t=' num2str(t) '\n']);
    %
    %           rev = 1;
    %       end
    
    
    %%
    
    n = n+1;
    
    % for methods that need a velocity field at n-1 the first time step
    % (e.g. AB-CN, oneleg beta) use ERK or IRK 

    if ((method_temp==2 || method_temp==5) && n<=method_startup_no)
        fprintf(fcw,['starting up with method ' num2str(method_startup) '\n']);
        method      = method_startup;        
    else
        method      = method_temp;
    end
    
    % perform a single time step with the time integration method
    if (method==2)
        [V,p,conv] = time_AB_CN(Vn,pn,conv_old,tn,dt,options);
        conv_old = conv;
    elseif (method==5)
        [V,p] = time_oneleg(Vn,pn,V_old,p_old,tn,dt,options);
    elseif (method==20)
        [V,p] = time_ERK(Vn,pn,tn,dt,options);
    elseif (method==21)
        [V,p,nonlinear_its(n)] = time_IRK(Vn,pn,tn,dt,options);
    else
        error('time integration method unknown');
    end
    
    
    % the velocities and pressure that are just computed are at
    % the new time level t+dt:
    t = tn + dt;
    time(n) = t;
                
    % update old solution (used for multistep methods)
    V_old = Vn;
    p_old = pn;    
    
    % update solution
    Vn = V;
    pn = p;
    tn = t;        
    
    % check residuals, conservation, set timestep, write output files
    process_iteration;
    
    
end