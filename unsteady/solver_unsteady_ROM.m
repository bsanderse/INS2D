%  Main solver file for unsteady calculations with reduced order model


%% load snapshot data

disp('loading data and making SVD...');
snapshots = load(snapshot_data,'uh_total','vh_total','dt','t_end','Re','k','umom','vmom');
% snapshots.U = [snapshots.uh_total; snapshots.vh_total];

% dt that was used for creating the snapshot matrix:
dt_snapshots = snapshots.dt;
% velocity field as snapshot matrix:
V_total      = [snapshots.uh_total';snapshots.vh_total'];

% check input dimensions
Nspace  = size(V_total,1); % total number of unknowns (Nu+Nv) of the original model
Nu      = options.grid.Nu;
Nv      = options.grid.Nv;
if (Nspace ~= Nu+Nv)
    error('The dimension of the snapshot matrix does not match the input dimensions in the parameter file');
end
 
if (snapshots.Re ~= Re)
    error('Reynolds numbers of snapshot data and current simulation do not match');
end


%% construct economic SVD
% note uh_total is stored as a Nt*Nu matrix instead of Nu*Nt which we use for
% the SVD

% subtract non-homogeneous BC contribution
if (options.rom.rom_bc == 1)
    f       = options.discretization.yM;
    dp      = pressure_poisson(f,t,options);
    Vbc     = - options.grid.Om_inv.*(options.discretization.G*dp);
    V_total = V_total - Vbc; % this velocity field satisfies M*V_total = 0
else
    Vbc = zeros(Nu+Nv,1);
end
options.rom.Vbc = Vbc;
    
    

% sample dt is multiple of snapshot dt:
if (rem(dt_sample,dt_snapshots) == 0)
    Nskip = dt_sample/dt_snapshots;
    % check if t_sample is multiple of dt_sample
    if (rem(t_sample,dt_sample) == 0)
        Nsnapshots    = t_sample / dt_snapshots; %size(V_total,2);
        snapshot_indx = 1:Nskip:Nsnapshots;
    else
        error('sample dt is not an integer multiple of sample time');
    end
else
    error('sample dt is not an integer multiple of snapshot dt');
end


% select snapshots
V_svd = V_total(:,snapshot_indx);

% enforce momentum conservation (works for periodic domains)
if (options.rom.mom_cons == 1)

    e_u = zeros(Nspace,1);
    e_v = zeros(Nspace,1);
    e_u(1:Nu)     = 1;
    e_v(Nu+1:end) = 1;
    e = [e_u e_v];
    e = e / norm(e);
    
    % 1) take SVD of (I-yy')*V_svd
    Vmod = V_svd - e*(e'*V_svd);
    [W,S,Z] = svd(Vmod,'econ');
    % 2) add y
    W = [e W];

%     disp('error in representing vector y before truncating:');
%     norm(Wext*Wext'*y - y,'inf')

    % 3) truncate:
%     Wextk = Vext(:,1:k);
%     Wmodk = Wmod(:,1:k);
%     Smodk = Smod(1:k,1:k);
%     disp('error in representing vector y after truncating:');
%     norm(Vextk*Vextk'*y - y,'inf')
% 
%     Umodk = Vextk*Smodk*Wmodk';
else
    % perform SVD
    [W,S,Z] = svd(V_svd,'econ');

end

% take first M columns of W as a reduced basis
% maximum:
% M = size(Wu,2);
% reduction:
% M = floor(Nspace/100);
% M = 16;
% options.rom.M = M;
% (better is to look at the decay of the singular values in S)
B  = W(:,1:M);
% Bu = Wu(:,1:M);
% Bv = Wv(:,1:M);
Bu = B(1:Nu,:);
Bv = B(Nu+1:end,:);
options.rom.B = B;
options.rom.Bu = Bu;
options.rom.Bv = Bv;
% options.rom.BuT = BuT;
% options.rom.BvT = BvT;
toc

% relative information content:
Sigma = diag(S);
RIC  = sum(Sigma(1:M).^2)/sum(Sigma.^2);
disp(['relative energy captured by SVD = ' num2str(RIC)]);
figure
semilogy(Sigma/Sigma(1),'s');
% or alternatively
% semilogy(Sigma.^2/sum(Sigma.^2),'s');

disp('starting time-stepping...');

%% precompute matrices
if (options.rom.precompute_diffusion == 1)
    Nu = options.grid.Nu;
    Nv = options.grid.Nv;
    options.rom.Diff  = Bu'*spdiags(options.grid.Omu_inv,0,Nu,Nu)*options.discretization.Diffu*Bu + ...
                        Bv'*spdiags(options.grid.Omv_inv,0,Nv,Nv)*options.discretization.Diffv*Bv;
    options.rom.yDiff = Bu'*spdiags(options.grid.Omu_inv,0,Nu,Nu)*options.discretization.yDiffu + ...
                        Bv'*spdiags(options.grid.Omv_inv,0,Nv,Nv)*options.discretization.yDiffv;
end

%% initialize reduced order solution
V  = V - Vbc; % subtract boundary condition contribution

R  = B'*V;

% map back to velocity space to get statistics of initial velocity field
% note that V will not be equal to the specified initial field, because
% B*B' does not equal identity in general
V  = B*R + Vbc; % add boundary condition contribution

[maxdiv(1), umom(1), vmom(1), k(1)] = check_conservation(V,t,options);


%% reduced order solution

% test implementation as follows:
% uh = rand(options.grid.Nu,1);
% vh = rand(options.grid.Nv,1);
% 
% ru = Bu'*uh;
% rv = Bv'*vh;
% 
% [convu_ROM,convv_ROM] = convectionROM(ru,rv,Bu,Bv,t,options,0);
% % this should equal using the original code but with B*r as input:
% [convu,convv] = convection(Bu*ru,Bv*rv,t,options,0);
% error_convu = Bu'*convu - convu_ROM;
% error_convv = Bv'*convv - convv_ROM;
% plot(error_convv)

%% load restart file if necessary
if (options.output.save_results==1)
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

%% for multistep methods or methods that need extrapolation of previous time steps
% if (method==2) % the only multistep method considered sofar
%     cu     = uh;
%     cv     = vh;
%     convection;
%     Cu_old = du2dx + duvdy;
%     Cv_old = duvdx + dv2dy;
% end
% 
% % for methods that need u^(n-1)
% uh_old = uh;
% vh_old = vh;
% 
% % for methods that need extrapolation of convective terms
% if (method == 62 || method == 92 || method==142 || method==172 || method==182 || method==192)
%     V_ep      = zeros(Nu+Nv,method_startup_no);
%     V_ep(:,1) = V;
% end


dtn    = dt;

eps    = 1e-12;

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
    %         Cu = Cux*spdiags(Iu_ux*uh+yIu_ux,0,N1,N1)*Au_ux + ...
    %              Cuy*spdiags(Iv_uy*vh+yIv_uy,0,N2,N2)*Au_uy;
    %         Cv = Cvx*spdiags(Iu_vx*uh+yIu_vx,0,N3,N3)*Av_vx + ...
    %              Cvy*spdiags(Iv_vy*vh+yIv_vy,0,N4,N4)*Av_vy;
    %
    %         Eu = eig(full(spdiags(Omu_inv,0,Nu,Nu)*Cu));
    %         Ev = eig(full(spdiags(Omv_inv,0,Nv,Nv)*Cv));
    %         max_eig(n) = max(abs([Eu;Ev]))*dt;
    %         dtn = dt;
    %         set_timestep;
    %         dt  = dtn;
    %         max_eig_G(n) = labda_conv*dt;
    %%
    
    n = n+1;
    
    % for methods that need a velocity field at n-1 the first time step
    % use RK4 or FC2: one-leg, AB, CN, FC1, IM
    % CN needs start-up for extrapolated Picard linearization
%     if ((method_temp==2 || method_temp==4 || method_temp==5 || ...
%             method_temp==71 || method_temp==62  || ...
%             method_temp==92 || method_temp==142 || method_temp==172  || method==182 || method==192)...
%             && n<=method_startup_no)
%         fprintf(fcw,['starting up with ' num2str(method_startup) '\n']);
%         method      = method_startup;
%         
%     else
%         method      = method_temp;
%     end
    
    
    % perform one time step with the time integration method
    if (method == 20)
        time_ERK_ROM;          
    elseif (method == 21)
        time_IRK_ROM;
    end

    
%     select_time_method;
    
    
    % the velocities and pressure that are just computed are at
    % the new time level t+dt:
    t = t + dt;
    time(n) = t;
    
    % check residuals, conservation, write output files
    process_iteration;
    
    
end