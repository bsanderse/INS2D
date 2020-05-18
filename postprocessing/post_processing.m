%% post-processing

file_name = [case_name '_pp.m'];
full_name = [folder_cases '/' case_name '/' file_name];

if (exist(full_name,'file'))
    
    run(full_name);
    
else
    
    disp(['postprocessing file ' file_name ' not available']);
    
end

%% additional Reduced-Order Model postprocessing
if (options.rom.rom == 1)
    % check if ROM simulation dt is same as FOM dt, or an integer multiple of
    % it
    if (rem(dt,dt_snapshots) == 0)
        skip = dt/dt_snapshots;
        % final time should be smaller than FOM time
        if (t_end<=snapshots.t_end)
            snapshot_end = ceil(t_end/dt_snapshots);
            snapshot_indx = 1:skip:(snapshot_end+1);
            t_vec = t_start:dt:t_end;
            
            % if velocity fields have been stored, we can compute errors
            if (options.rom.process_iteration_FOM==1)
                
                if (options.output.save_unsteady == 1)
                    % we have the velocity fields, so we can compute error wrt
                    % FOM
                    % uh_total is of size Nt*Nu, V_total size (Nu+Nv)*Nt
                    V_total = [uh_total vh_total]';
                    snapshots_V_total = [snapshots.uh_total(snapshot_indx,:) snapshots.vh_total(snapshot_indx,:)]';
                    error_V = V_total - snapshots_V_total;

                    % inf-norm
                    error_V_inf = max(abs(error_V),[],1);
                    
                    % 2-norm of error
                    % note that 2-norm of velocity-field is simply sqrt(2*k),
                    % with k the kinetic energy = 0.5*V'*Om*V
                    % NOTE! the (finite volume)-weighted 2-norm is consistent with the
                    % Frobenius norm of the optimization problem as solved
                    % by the SVD                    
%                     Om_total = sum(sum(options.grid.Om));
                    
                    % choose reference velocity field (can be time
                    % dependent)
                    % V_ref = 1:
                    V_ref    = ones(size(snapshots_V_total));
                    % V_ref = snapshots:
%                     V_ref    = snapshots_V_total;
                    
                    V_2_ref  = weightedL2norm(V_ref,options.grid.Om); 
                    error_V_2 = weightedL2norm(error_V,options.grid.Om)./V_2_ref;
                    
                    
                    % best possible approximation given the projection:
                    % note that the norm should be consistent with the
                    % optimization problem used in the SVD
                    V_best = getFOM_velocity(getROM_velocity(snapshots_V_total,0,options),0,options);
                    error_V_best = V_best - snapshots_V_total;
                    error_V_best_2 = weightedL2norm(error_V_best,options.grid.Om)./V_2_ref;
                    
                    figure(101)
%                     plot(t_vec,error_V_inf);
%                     hold on
                    % skip i=1, as error_v_2_norm is zero for i=1
                    plot(t_vec,error_V_2); %(2:end)./error_V_2_norm(2:end));                    
                    hold on
                    plot(t_vec,error_V_best_2);
                    set(gca,'Yscale','log');
                    legend('L_2 error in ROM velocity','Best approximation (projection FOM)')
%                     legend('L_{inf} error in ROM velocity','L_2 error in ROM velocity','Best approximation (projection FOM)')

                    
                    if (options.rom.pressure_recovery == 1)
                        % correct spatial mean of both to be zero
                        % p_total is of size Nt*Np, change to Np*Nt
                        p_total = p_total';
                        mean_ROM = mean(p_total,1);
                        snapshots_p_total = snapshots.p_total(snapshot_indx,:)';
                        mean_FOM = mean(snapshots_p_total,1);
                        
                        error_p = (p_total - mean_ROM) - (snapshots_p_total - mean_FOM);

                        % inf-norm
                        error_p_inf = max(abs(error_p),[],1);
                        
                        % 2-norm of error    
                        % choose reference pressure field (can be time
                        % dependent)                        
                        p_ref     = 0.25*ones(size(snapshots_p_total));
                        % p_ref = snapshots:
%                         p_ref    = snapshots_p_total - mean_FOM;

                        p_2_ref   = weightedL2norm(p_ref,options.grid.Omp); 
                        error_p_2 = weightedL2norm(error_p,options.grid.Omp)./p_2_ref;

                        % best possible approximation given the projection:                        
                        p_best = getFOM_pressure(getROM_pressure(snapshots_p_total,0,options),0,options);
                        mean_ROM_best = mean(p_best,1);
                        
                        error_p_best = (p_best - mean_ROM_best) - (snapshots_p_total - mean_FOM);
                        error_p_best_2 = weightedL2norm(error_p_best,options.grid.Omp)./p_2_ref;
                        
                        figure(102)
                        plot(t_vec,error_p_2);
                        hold on
%                         plot(t_vec,error_p_inf);
                        plot(t_vec,error_p_best_2);
                        set(gca,'Yscale','log');
                        legend('L_{2} error in ROM pressure','Projection FOM pressure')

                    end
                    

                end
                
                
                figure(103)
%                 semilogy(t_vec,abs(k - snapshots.k(snapshot_indx))/snapshots.k(1));
%                 semilogy(t_vec,abs(k - snapshots.k(1))/snapshots.k(1));
                semilogy(t_vec,abs(k-k(1))/k(1));
                hold on
%                 set(gca,'Yscale','log')
                ylabel('energy error');

%                 legend('(K_{ROM}(t)-K_{FOM}(t))/K_{FOM}(0)','(K_{ROM}(t)-K_{FOM}(0))/K_{FOM}(0)','(K_{ROM}(t)-K_{ROM}(0))/K_{ROM}(0)')
%                 title('error in kinetic energy ROM');
                
                figure(104)
%                 umom0 = snapshots.umom(1);
                umom0 = umom(1);
                semilogy(t_vec,abs(umom-umom0)/umom0)
                hold on
                ylabel('momentum error');
                
                figure(105)
                plot(t_vec,maxdiv);
                hold on
                plot(t_vec,snapshots.maxdiv(snapshot_indx));
                xlabel('t')
                ylabel('maximum divergence of velocity field');
                grid
                
                
            end
            %         hold on
            %         plot(t_start:dt:t_end,error_v);
            %         legend('error in u','error in v');
            
        end
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
