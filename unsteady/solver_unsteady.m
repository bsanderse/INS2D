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
    [convu,convv] = convection(V,V,t,options,0);
    rhs_terms.conv = [convu;convv];
    
    switch options.case.boussinesq
        case 'temp'
            rhs_terms.convT = convection_temperature(T,V,t,options,0);
            switch options.temp.incl_dissipation
                case 1
                    rhs_terms.Phi = dissipation(V,t,options,0);
            end
            
    end
    
end

% for methods that need u^(n-1)
if (method == 5)
    V_old  = V;
    p_old  = p;
    T_old  = T;
end

% set current velocity and pressure
Vn = V;
pn = p;
Tn = T; % temperature
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
time_start = toc

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
    switch method
        case 2
            [V,p,T,rhs_terms] = time_AB_CN(Vn,pn,Tn,rhs_terms,tn,dt,options);
            %         conv_old = conv;
            %         convT_old = convT;
        case 5
            [V,p] = time_oneleg(Vn,pn,V_old,p_old,tn,dt,options);
        case 20
            [V,p,T] = time_ERK(Vn,pn,Tn,tn,dt,options);
        case 21
            [V,p,nonlinear_its(n)] = time_IRK(Vn,pn,tn,dt,options);
        case 22
            [V,p,T,nonlinear_its(n)] = time_IM_Boussinesq(Vn,pn,Tn,tn,dt,options);
        case 23
            [V,p,T,nonlinear_its(n)] = time_IM_Boussinesq_split(Vn,pn,Tn,tn,dt,options);
        otherwise
            error('time integration method unknown');
    end
    
    
    % the velocities and pressure that are just computed are at
    % the new time level t+dt:
    t = tn + dt;
    time(n) = t;
    
    % check residuals, conservation, set timestep, write output files
    process_iteration;
    
    % update solution
    V_old = Vn;
    p_old = pn;
    T_old = Tn;
    Vn = V;
    pn = p;
    Tn = T;
    tn = t;
    
    
    
end
disp('finished time-stepping...');
time_loop = toc-time_start