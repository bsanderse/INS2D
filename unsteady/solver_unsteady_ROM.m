%  Main solver file for unsteady calculations with reduced order model


solver_unsteady_ROM_basis_construction;

%% precompute ROM operators by calling operator_rom
% results are stored in options structure
disp('precomputing ROM operators...');
precompute_start = toc;
options = operator_rom(options);
precompute_end(j) = toc-precompute_start


%% initialize reduced order solution
% we expand the part of the solution vector that is div-free in terms of
% B*R
% V = B*R + Vbc
V_orig = V;
% get the coefficients of the ROM
R = getROM_velocity(V,t,options);

if (options.rom.process_iteration_FOM == 1)
    % map back to velocity space to get statistics of initial velocity field
    % note that V will not be equal to the specified initial field, because
    % B*B' does not equal identity in general
    V  = getFOM_velocity(R,t,options);
    
    [maxdiv(1), umom(1), vmom(1), k(1)] = check_conservation(V,t,options);
    
    
    % overwrite the arrays with total solutions
    if (save_unsteady == 1)
        uh_total(n,:) = V(1:options.grid.Nu);
        vh_total(n,:) = V(options.grid.Nu+1:end);
    end
end


if (options.rom.div_free == 0)
    % get initial ROM pressure by projecting initial presssure
    q = getROM_pressure(p,t,options);
    
elseif (options.rom.div_free == 1)
    if (options.rom.pressure_recovery == 1)
        % get initial pressure that corresponds to the ROM velocity field
        q = pressure_additional_solve_ROM(R,t,options);
        p = getFOM_pressure(q,t,options);

        % overwrite the arrays with total solutions
        if (save_unsteady == 1)
            p_total(n,:)  = p;
        end
    else
        % q is not used, set to arbitrary value
        q = 0;
    end
end


%% reduced order solution tests

% % test the offline implementation as follows:
% Rtest = rand(options.rom.M,1);
% Vtest = getFOM_velocity(Rtest,t,options);
%
% % with precomputing:
% conv_ROM_pre = convectionROM(Rtest,t,options,0);
% diff_ROM_pre = diffusionROM(Rtest,t,options,0);
%
% % without precomputing:
% [convu, convv, dconvu, dconvv] = convection(Vtest,Vtest,t,options,0);
% conv_ROM  = B'*(Om_inv.*[convu;convv]);
% [d2u,d2v,dDiffu,dDiffv] = diffusion(Vtest,t,options,0);
% diff_ROM  = B'*(Om_inv.*[d2u;d2v]);
%
% % compute error between the two versions
% error_conv_ROM = conv_ROM_pre - conv_ROM;
% error_diff_ROM = diff_ROM_pre - diff_ROM;
% plot(error_conv_ROM)
% hold on
% plot(error_diff_ROM);

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


%% start time stepping

dtn    = dt;

eps    = 1e-12;

disp('starting time-stepping...');

time_start = toc

while(n<=nt)
    
    %% dynamic timestepping:
    % set_timestep;
    
    %%
    
    n = n+1;
    
    % time reversal (used in inviscid shear-layer testcase)
    %     if (abs(t-t_end/2)<1e-12)
    %         dt = -dt;
    %     end
    
    % perform one time step with the time integration method
    if (method == 20)
        time_ERK_ROM;
    elseif (method == 21)
        time_IRK_ROM;
    else
        error('time integration method unknown');
    end
    
    
    % the velocities and pressure that are just computed are at
    % the new time level t+dt:
    t = t + dt;
    time(n) = t;
    
    % check residuals, conservation, write output files
    % this requires to go back to FOM level
    if (options.rom.process_iteration_FOM == 1)
        % map back from reduced space to full order model space
        % this is used for postprocessing purposes, e.g. evaluating the divergence
        % of the velocity field
        V = getFOM_velocity(R,t,options);
        
        process_iteration;
    end
    
    
end
disp('finished time-stepping...');
time_loop(j) = toc-time_start


V  = getFOM_velocity(R,t,options);