%  Main solver file for unsteady calculations with reduced order model

switch options.rom.rom_type
    
    case 'POD'
        %% POD
        % load snapshot data
        % assume that for parametric studies (e.g. changing number of modes), the
        % FOM data file does not change, so we only load it once
        if (j==1)
            disp(['loading datafile...: ' snapshot_data]);
            snapshots = load(snapshot_data,'uh_total','vh_total','p_total','T_total','dt','t_end','Re','k','umom','vmom','maxdiv','Vbc');

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
        
        switch options.case.boussinesq   
            case 'temp'
                options.rom.BT = getTemperatureBasisPOD(snapshots,snapshot_sample_index,options); % temperature rom basis
        end
    otherwise
        error ('wrong ROM type')
        
end
