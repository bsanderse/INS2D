%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% create files and directories

if (restart.load==0 && options.output.save_results==1)
    file         = [project,'_',num2str(Re,'%10.3e'),'_',num2str(L_x),'x',num2str(L_y),...
        '_',num2str(Nx),'x',num2str(Ny)];
    
    % create results directory if not present
    if (~isdir(path_results))
        mkdir(path_results);
    end
    
    path_results = [path_results '/' file];
    
    i_file = 0;
    while (exist(path_results,'dir')==7)
        path_results = ['results/' file '_' num2str(i_file)];
        i_file = i_file+1;
    end
    mkdir(path_results);
    
    % copy inputfiles in result folder so that restart is possible and settings
    % can be reproduced
    %     mkdir([path_results '/inputfiles']);
    %     copyfile('inputfiles/*.m ' path_results '/inputfiles']);
    
end

% newline and tab
nl         = sprintf('\n');
tab        = sprintf('\t');

% filenames
if (options.output.save_results == 1)
    
    file_pres  = [path_results,'/iterations_pressure.txt'];
    file_conv  = [path_results,'/convergence.txt']; % convergence / timestep information
    file_cw    = [path_results,'/cw_output.txt']; % command window output
    
    if (restart.load == 0)
        if (steady==0)
            % file with information pressure solve
            fpres      = fopen(file_pres,'w+');
            fwrite(fpres,[project ' ' num2str(Nx) ' ' num2str(Ny) ' ' num2str(Re,'%10.3e') nl]);
            fwrite(fpres,['convergence information pressure solve' nl]);
            fprintf(fpres,'n          iter.        norm             cpu-time\n');
        else
            fpres      = 0;
        end
        
        % file with convergence information
        fconv      = fopen(file_conv,'w+');
        fwrite(fconv,[project ' ' num2str(Nx) ' ' num2str(Ny) ' ' num2str(Re,'%10.3e') nl]);
        fwrite(fconv,['residual and conservation information' nl]);
    else
        % open such that data is appended
        fpres      = fopen(file_pres,'a+');
        fconv      = fopen(file_conv,'a+');
    end
    
    
    if (cw_output == 0)
        fcw        = fopen(file_cw,'a+');
    else
        fcw        = 1;
    end
    

    
else
    fcw   = 1;
    fpres = 0;
    fconv = 0;
%     options.output.fcw = fcw;
end
options.output.fpres = fpres;
options.output.fconv = fconv;
options.output.fcw   = fcw;

% save workspace to a .mat file at the end of simulation
if (options.output.save_file == 1)
    file_mat   = [path_results,'/matlab_data.mat'];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%