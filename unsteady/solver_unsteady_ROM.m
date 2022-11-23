%  Main solver file for unsteady calculations with reduced order model


switch options.rom.rom_type
    
    case 'POD'
        %% POD
        % load snapshot data
        % assume that for parametric studies (e.g. changing number of modes), the
        % FOM data file does not change, so we only load it once
        if (j==1)
            disp(['loading datafile...: ' snapshot_data]);
            snapshots = load(snapshot_data,'uh_total','vh_total','p_total','dt','t_end','Re','k','umom','vmom','maxdiv','Vbc');

            % dt that was used for creating the snapshot matrix:
            dt_snapshots = snapshots.dt;
            options.rom.dt_snapshots = dt_snapshots;
            
            if (snapshots.Re ~= Re)
                error('Reynolds numbers of snapshot data and current simulation do not match');
            end
            
            % find indices of snapshot matrix that are needed
            snapshot_sample_index = getSampleIndex(options.rom);
           
        end
        
        % construct velocity basis        
        [options.rom.B, options.rom.div_free, options.rom.Vbc, options.rom.yM] = getVelocityBasisPOD(snapshots,snapshot_sample_index,options);
        % construct pressure basis
        if (options.rom.pressure_recovery == 1 || options.rom.div_free == 0)
            [options.rom.Bp,options.rom.Mp] = getPressureBasisPOD(snapshots,snapshot_sample_index,options);
        end
        
    case 'Fourier'

        Vbc    = zeros(options.grid.NV,1);
        options.rom.Vbc  = Vbc;
        options.rom.yM   = 0;
        options.rom.div_free = 0;
        [options.rom.B, options.rom.Bp, options.rom.M] = getVelocityBasisFourier(options);

        
    case 'FDG'
        
        Vbc    = zeros(options.grid.NV,1);
        options.rom.Vbc  = Vbc;
        options.rom.yM   = 0;
        options.rom.div_free = 1;
        [options.rom.B, options.rom.B_inv, options.rom.M] = getVelocityBasisFDG(options);
        
    case 'FDG-Fourier'
        
        
        Vbc    = zeros(options.grid.NV,1);
        options.rom.Vbc  = Vbc;
        options.rom.yM   = 0;
        options.rom.div_free = 1;
        [options.rom.B, options.rom.B_inv, options.rom.M] = getVelocityBasisFDGFourier(options);
                
    otherwise
        error ('wrong ROM type')
        
end

%% precompute ROM operators by calling operator_rom
% results are stored in options structure
disp('precomputing ROM operators...');
precompute_start = toc;
options = operator_rom(options);
precompute_end(j) = toc-precompute_start


%% initialize reduced order solution
[R,q] = initializeROM(V,p,t,options);


% map back to FOM space to get initial properties
if (options.rom.process_iteration_FOM == 1)
    % map back to velocity space to get statistics of initial velocity field
    % note that V will not be equal to the specified initial field, because
    % B*B' does not equal identity in general
    V  = getFOM_velocity(R,t,options);
    if (options.rom.pressure_recovery == 1)
        p = getFOM_pressure(q,t,options);
    end
    [maxdiv(1), umom(1), vmom(1), k(1)] = check_conservation(V,T,t,options);
    
    
    % overwrite the arrays with total solutions
    if (save_unsteady == 1)
        uh_total(n,:) = V(options.grid.indu);
        vh_total(n,:) = V(options.grid.indv);
        if (options.rom.pressure_recovery == 1)
            p_total(n,:)  = p;
        end
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

% set current velocity and pressure
Rn = R;
qn = q;
tn = t;

eps    = 1e-12;

disp('starting time-stepping...');

time_start = toc

while(n<=nt)
    
    %% dynamic timestepping:
    % set_timestep;
    
    %% perform one time step with the time integration method
    
    n = n+1;
       
    if (method == 20)
         [R,q] = time_ERK_ROM(Rn,qn,tn,dt,options);
    elseif (method == 21)
        [R,q,nonlinear_its(n)] = time_IRK_ROM(Rn,qn,tn,dt,options);
    else
        error('time integration method unknown');
    end    
    
    % the velocities and pressure that are just computed are at
    % the new time level t+dt:
    t       = tn + dt;
    % store current time in vector
    time(n) = t;
    
    % update old solution
    Rn = R;
    qn = q;
    tn = t;   
    
    
    % check residuals, conservation, write output files
    % this requires to go back to FOM level
    if (options.rom.process_iteration_FOM == 1)
        % map back from reduced space to full order model space
        % this is used for postprocessing purposes, e.g. evaluating the divergence
        % of the velocity field
        V = getFOM_velocity(R,t,options);
        
        if (options.rom.pressure_recovery == 1)
            p = getFOM_pressure(q,t,options);
        end
        process_iteration;
    end
    

end
disp('finished time-stepping...');
time_loop(j) = toc-time_start

% get FOM velocity and pressure for postprocessing purposes
V  = getFOM_velocity(R,t,options);
if (options.rom.pressure_recovery == 1)
    p = getFOM_pressure(q,t,options);
end