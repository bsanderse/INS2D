if (exist('nonlinear_its'))
nonlinear_its_avg(j,jj) = mean(nonlinear_its(2:end));
disp(['average number of nonlinear iterations: ' num2str(nonlinear_its_avg(j,jj))]);
end

% if (jj==1)
%     
%     ufine = uh;
%     vfine = vh;
%     pfine = p;
%     
% end

% save(['results/vanKan/implicit/vanKan_RadauIIAB_N20_t1_dt' num2str(dt) '.mat']);
% keyboard;
% if (jj>1 || length(dt_list)==1)
    load(['results/vanKan/implicit/vanKan_ERK4_N' num2str(Nx) '_t1.mat'],'ufine','vfine','pfine');
        
    

    % temporal error obtained by subtracting spatial error (obtained with
    % ufine)
    u_error = uh-ufine(:);
    v_error = vh-vfine(:);
    
    p_error = p - pfine(:);
    
    
    u_error_i(j,jj) = max(abs(u_error))
    u_error_2(j,jj) = sqrt(sum( u_error.^2)/Nu)
    v_error_i(j,jj) = max(abs(v_error))
    v_error_2(j,jj) = sqrt(sum( v_error.^2)/Nv)
    p_error_i(j,jj) = max(abs(p_error))
    p_error_2(j,jj) = sqrt(sum( p_error.^2)/Np)
       
%     uspace_error(j) = max(abs(ufine_error));
%     vspace_error(j) = max(abs(vfine_error));
%     pspace_error(j) = max(abs(pfine_error));
% end

if (jj==length(dt_list) && j==length(mesh_list))
    
%     % spatial order
%     disp('spatial order:')
%     get_slope(1./mesh_list,abs(uspace_error))
%     get_slope(1./mesh_list,abs(vspace_error))
%     get_slope(1./mesh_list,abs(pspace_error))
%     figure(2)
%     loglog(1./mesh_list,abs(uspace_error),'bx-');
%     hold on
%     loglog(1./mesh_list,abs(vspace_error),'rx-');
%     loglog(1./mesh_list,abs(pspace_error),'kx-');
    
    colors1 = {'kx-','rx-','bx-','gx-','mx-','yx-','cx-'};
    colors2 = {'ko--','ro--','bo--','go--','mo--','yo--','co--'};
    
    for k = 1:length(mesh_list)
        figure(1)
        loglog(dt_list(2:end),abs(u_error_i(k,2:end)),char(colors1(k)))
        hold on
        loglog(dt_list(2:end),abs(u_error_2(k,2:end)),char(colors2(k)))
        figure(2)
        loglog(dt_list(2:end),abs(p_error_i(k,2:end)),char(colors1(k)))
        hold on
        loglog(dt_list(2:end),abs(p_error_2(k,2:end)),char(colors2(k)))        
        
        % temporal order for each grid, using errors at finest dt
        get_slope(dt_list(2:end),abs(u_error_i(k,2:end)),1,2)
        get_slope(dt_list(2:end),abs(u_error_2(k,2:end)),1,2)
        get_slope(dt_list(2:end),abs(p_error_i(k,2:end)),1,2)
        get_slope(dt_list(2:end),abs(p_error_2(k,2:end)),1,2)
    
    end
    
    legend('Linf N=5','L2 N=5','Linf N=10','L2 N=10','Linf N=20','L2 N=20',...
           'Linf N=40','L2 N=40');
    
    % average order reduction, assuming mesh is halved each time; take
    % results at third finest time step
    disp('order reduction: ');
    log(abs(u_error_i(2:end,3))./abs(u_error_i(1:end-1,3)))/log(2)
    log(abs(u_error_2(2:end,3))./abs(u_error_2(1:end-1,3)))/log(2)
    
   
end