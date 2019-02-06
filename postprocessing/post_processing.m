%% post-processing

file_name = [case_name '_pp.m'];
full_name = [folder_cases '/' case_name '/' file_name];

if (exist(full_name,'file'))
    
    run(full_name);

else
    
    disp(['postprocessing file ' file_name ' not available']);
    
end

%% roms
if (options.rom.rom == 1)
    % check if ROM simulation dt is same as FOM dt, or an integer multiple of
    % it
    if (rem(dt,dt_snapshots) == 0)
        skip = dt/dt_snapshots;
        snapshot_indx = 1:skip:size(snapshots.uh_total,1);
        % uh_total is of size Nt*Nu
        error_u = max(abs(uh_total - snapshots.uh_total(snapshot_indx,:)),[],2);
        error_v = max(abs(vh_total - snapshots.vh_total(snapshot_indx,:)),[],2);
        figure
        plot(t_start:dt:t_end,error_u);
        hold on
        plot(t_start:dt:t_end,error_v);
        legend('error in u','error in v');
        title('error in ROM velocity components');
    end
    
end

%% standard plots
% velocity;
% 
% pressure;
% 
% vorticity;
% 
% streamfunction;

%% Poiseuille
% u = reshape(uh,Nux_in,Nuy_in);
% u_exact = -6*yp.*(yp-L_y);
% 
% e2(j) = sqrt(sum((u(1,:)'-u_exact).^2)/Npy)
% einf(j) = max(abs(u(1,:)'-u_exact))
% p_add_solve=1;
% pressure_additional_solve;
% run('results/poiseuille_blowing/pressure_study/error_finalp.m');

%% Shih LDC
% pp_LDC_Shih;
% pp_LDC_unsteady;

%% LDC_unsteady
% k_list(j) = k(end);
% pp_LDC;
% p_add_solve=1;
% pressure_additional_solve;
% run('results/LDC/pressure_study/error_finalp.m');

%% BFS
% pp_BFS;

%% AD
% pp_AD;
% velocity
% pressure

%     uh = uh+cos(t); % from moving ref to stationary ref.
% 
%     u_curr   = reshape(uh,Nux_in,Nuy_in);
%     xin_curr = xin;
%     yp_curr  = yp;
%     load('../v2.4/results/AD_2D_translating/periodic_domain/N512_movingref_new.mat','uh','xin','yp');
% %     
%     ufine_int = interp2(xin,yp',reshape(uh,length(xin),length(yp)),xin_curr,yp_curr','spline');
%     error(j) = max(abs(u_curr(:)-ufine_int(:)))
%     error2(j) = sqrt(sum((u_curr(:)-ufine_int(:)).^2)/Nu)
    

% velocity
% line = {'rx-','bo-','ms-','kd-','g-','y-'};
% color = char(line(j));
% wake_pos = find(xin==1);
% % % disp(['position wake: ' num2str(xin(wake_pos))])
% % % 
% figure(1)
% % uh = uh+cos(t); % from moving ref to stationary ref.
% 
% u = reshape(uh,Nux_in,Nuy_in);
% velocity
% pressure
% pause
% u_wake = u(wake_pos,:);
% if (j==1)
%     load('results/AD_2D_translating/periodic_domain/N512_movingref_wake_Re100_order2.mat','ufine','yfine');
%     ufine = u_wake; 
%     yfine = yp;
%     save(['results/AD_2D_translating/periodic_domain/N' num2str(Nx) '_movingref_wake_Re100_order4.mat'],...
%         'uh','vh','p','Nx','Ny','xp','yp','xin','yin','t','dt','Re','ufine','yfine');   
% else

% if (jj>1)
%     u_curr   = reshape(uh,Nux_in,Nuy_in);
%     xin_curr = xin;
%     yp_curr  = yp;
%     load('results/AD_2D_translating/periodic_domain/N256_movingref_new_Re10.mat','uh','xin','yp');
%     
%     ufine_int = interp2(xin,yp',reshape(uh,length(xin),length(yp)),xin_curr,yp_curr','spline');
%     error(jj) = max(abs(u_curr(:)-ufine_int(:)))
%     error2(jj) = sqrt(sum((u_curr(:)-ufine_int(:)).^2)/Nu)
%     
% end
%     ufine_int = interp1(yfine,ufine,yp,'spline');
%     error(j) = max(abs(ufine_int'-u_wake))    
%     error2(j) = sqrt(sum((ufine_int'-u_wake).^2)/length(u_wake))

% end
% 
% plot(yp,u(wake_pos,:),color)
% hold on


%% doublejet
% line = {'r-','b-','k-','m-','g-'};
% color = char(line(j));
% figure(1)
% plot(0:dt:t_end,(enstrophy-enstrophy(1))/enstrophy(1),color)
% hold on
% if (j==1)
% % enstrophy_fine = enstrophy(end);
% k_ex = k(end);
% end
% error_enstrophy(j) = enstrophy(end)-enstrophy_fine
% k_ex   = (2/15)*pi^2*(13*exp(15/2)+17*exp(-15/2))/(exp(15/2)+exp(-15/2))+(1/400)*pi^2;
% k21(j) = k2(end)-k_ex;
% k11(j) = k(end)-k_ex;
% vort_ex = -12.56636293; % exact initial value
% vort_ex = -12.536625993619506; % 2nd order 160x160
% vort_ex = -12.536374311199099; % 4th order 80x80
% error_vort(j) = omega_total3(1)-vort_ex
% pp_shearlayer;


%% moving reference frame
% t_plot = dt:dt:t_end;
% figure
% udiff = max(abs(2+cos(t_plot)-umin(2:end)),abs(2+cos(t_plot)-umax(2:end)));
% plot(t_plot,udiff,'b')

%% wake studies
% pp_wake;

%% poiseuille-like flow
% p_add_solve=1;
% pressure_additional_solve;
% run('results/periodic_channel/error_finalp.m')

%% energy and reversibility studies
% pp_reversibility;

%% Taylor vortex
% pp_Taylor;
% pp_Taylor_spacetime;

% p_add_solve=1;
% pressure_additional_solve;
% run('results/Taylor_Green/pressure_study/error_finalp.m');


%% van Kan
% pp_vanKan;
% p_add_solve=1;
% pressure_additional_solve;
% run('results/cornerflow_ben/pressure_study/error_finalp.m');

%% Kovasznay
% pp_Kovasznay;

%% wake meandering
% pp_meandering;

%% actuators
% figure(1)
% hold on
% plot([xp(floor(Nx/5)) xp(floor(Nx/5))],[yp(yrange(1)) yp(yrange(end))])
% plot([xp(2*floor(Nx/5)) xp(2*floor(Nx/5))],[yp(yrange(1)+floor(Ny/4)) yp(yrange(end)+floor(Ny/4))],'k-')
% plot([xp(2*floor(Nx/5)) xp(2*floor(Nx/5))],[yp(yrange(1)-floor(Ny/4)) yp(yrange(end)-floor(Ny/4))],'k-')
% plot([xp(3*floor(Nx/5)) xp(3*floor(Nx/5))],[yp(yrange(1)) yp(yrange(end))],'k-')
% plot([xp(4*floor(Nx/5)) xp(4*floor(Nx/5))],[yp(yrange(1)-floor(Ny/4)) yp(yrange(end)-floor(Ny/4))],'k-')
% plot([xp(4*floor(Nx/5)) xp(4*floor(Nx/5))],[yp(yrange(1)+floor(Ny/4)) yp(yrange(end)+floor(Ny/4))],'k-')
