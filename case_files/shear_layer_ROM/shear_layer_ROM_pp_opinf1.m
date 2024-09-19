if (options.rom.rom == 1 && strcmp(options.rom.rom_type,'POD'))
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
                    
                    figure(7101)
%                     plot(t_vec,error_V_inf);
%                     hold on
                    % skip i=1, as error_v_2_norm is zero for i=1
                    label = "M = " + num2str(M)+options.rom.opinf_type;
                    plot(t_vec,error_V_2,"displayname", "L_2 velocity error " + label); %(2:end)./error_V_2_norm(2:end));                    
                    hold on
                    plot(t_vec,error_V_best_2, "displayname", "best approx error "+ label);
                    set(gca,'Yscale','log');
                    % legend('L_2 error in ROM velocity','Best approximation (projection FOM)')
%                     legend('L_{inf} error in ROM velocity','L_2 error in ROM velocity','Best approximation (projection FOM)')
                    legend('show')
                                    

                end
                
                
                figure(7103)
%                 semilogy(t_vec,abs(k - snapshots.k(snapshot_indx))/snapshots.k(1));
%                 semilogy(t_vec,abs(k - snapshots.k(1))/snapshots.k(1));
                semilogy(t_vec,abs(k-k(1))/k(1),"displayname",label);
                hold on
%                 set(gca,'Yscale','log')
                ylabel('energy error');
                % legend("M=2, implicit","M=4, implicit","M=8, implicit","M=16, implicit","M=2, explicit","M=4, explicit","M=8, explicit","M=16, explicit")
                legend('show');

%                 legend('(K_{ROM}(t)-K_{FOM}(t))/K_{FOM}(0)','(K_{ROM}(t)-K_{FOM}(0))/K_{FOM}(0)','(K_{ROM}(t)-K_{ROM}(0))/K_{ROM}(0)')
%                 title('error in kinetic energy ROM');
            end
        end
    end
end