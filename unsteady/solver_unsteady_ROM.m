%  Main solver file for unsteady calculations with reduced order model

%% load snapshot data
snapshots = load(snapshot_data,'uh_total','vh_total');
% snapshots.U = [snapshots.uh_total; snapshots.vh_total];

%% construct economic SVD
% [W,S,Z] = svd(snapshots.U,'econ');
[Wu,Su,Zu] = svd(snapshots.uh_total','econ');
[Wv,Sv,Zv] = svd(snapshots.vh_total','econ');
%
N = size(snapshots.uh_total,2);
% take first m columns of W
m = floor(N/20);
% (better is to look at the decay of the singular values in S)
Bu = Wu(:,1:m);
Bv = Wv(:,1:m);

%% precompute matrices
options.discretization.Diffu  = Bu'*options.discretization.Diffu*Bu;
options.discretization.yDiffu = Bu'*options.discretization.yDiffu;
options.discretization.Diffv  = Bv'*options.discretization.Diffv*Bv;
options.discretization.yDiffv = Bv'*options.discretization.yDiffv;


%% reduced order solution

% test implementation as follows:
uh = rand(options.grid.Nu,1);
vh = rand(options.grid.Nv,1);

ru = Bu'*uh;
rv = Bv'*vh;

[convu_ROM,convv_ROM] = convectionROM(ru,rv,Bu,Bv,t,options,0);
% this should equal using the original code but with B*r as input:
[convu,convv] = convection(Bu*ru,Bv*rv,t,options,0);
error_convu = Bu'*convu - convu_ROM;
error_convv = Bv'*convv - convv_ROM;
plot(error_convv)

%% load restart file if necessary
% if (restart.load == 0 && options.output.save_results==1)
%     fprintf(fconv,'n            dt               t                res              maxdiv           umom             vmom             k\n');
%     fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n',...
%         n,dt,t,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
% end

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
    select_time_method;
    
    
    % the velocities and pressure that are just computed are at
    % the new time level t+dt:
    t = t + dt;
    time(n) = t;
    
    % check residuals, conservation, write output files
    process_iteration;
    
    
end