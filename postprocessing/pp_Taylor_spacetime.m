% %% coarse and fine simulation
% if (jj==1)
%     
%     ufine = uh;
%     vfine = vh;
%     pfine = p-mean(p(:));
%     
% else
% 
%     % temporal error obtained by subtracting spatial error (obtained with
%     % ufine)
%     u_error = uh-ufine;
%     v_error = vh-vfine;
%     
%     p       = p-mean(p(:));
%     p_error = p-pfine;
%     
%     u_error_i(j) = max(abs(u_error))
%     v_error_i(j) = max(abs(v_error))
%     u_error_2(j) = sqrt( sum(u_error.^2)/Nu )
%     v_error_2(j) = sqrt( sum(v_error.^2)/Nv )
% 
%     p_error_i(j) = max(abs(p_error))
%     p_error_2(j) = sqrt( sum(p_error.^2)/Np)    
%     
%     
% end

%% coarse and fine simulation
% if (j==1)
%     
%     ufine = uh;
%     vfine = vh;
%     pfine = p-mean(p(:));
%     
% else
% 
%     % temporal error obtained by subtracting spatial error (obtained with
%     % ufine)
%     u_error = uh-ufine;
%     v_error = vh-vfine;
%     
%     p       = p-mean(p(:));
%     p_error = p-pfine;
%     
%     u_error_i(j-1) = max(abs(u_error))
%     v_error_i(j-1) = max(abs(v_error))
%     u_error_2(j-1) = sqrt( sum(u_error.^2)/Nu )
%     v_error_2(j-1) = sqrt( sum(v_error.^2)/Nv )
% 
%     p_error_i(j-1) = max(abs(p_error))
%     p_error_2(j-1) = sqrt( sum(p_error.^2)/Np)    
%     
% end
% 
% if (j==length(dt_list))
%     % average slope over entire interval
%     beta = get_slope(dt_list(2:end),u_error_i,1,length(dt_list)-1)
%     % error constant
%     u_error_i./(dt_list(2:end).^beta)
%     % based on last simulation:
%     beta = get_slope(dt_list(2:end),u_error_i,1,2)
%     C2(jj)   = u_error_i(1)./(dt_list(1).^beta)
%     u_error_newi(jj) = u_error_i(3);
%     u_error_new2(jj) = u_error_2(3);
%     p_error_newi(jj) = p_error_i(3);
%     p_error_new2(jj) = p_error_2(3);
% 
%     
%     figure(1)
%     loglog(dt_list(2:end),u_error_i,'bx-')
%     hold on
% %     loglog(dt_list(2:end),u_error_2,'bx--')
% 
%     clear u_error_i;
% end

%% error wrt exact solution
 
t_factor1 = exp(-4*pi^2*t_end/Re);
t_factor2 = exp(-2*pi^2*t_end/Re);

% point value
u_exact = -sin(pi*xu).*cos(pi*yu)*t_factor2;
v_exact = cos(pi*xv).*sin(pi*yv)*t_factor2;

% volume average
% test = kron(y,ones(Nux_in,1));
% test = reshape(test,Nux_in,Ny+1);
% u_exact = 1/(hx(1)*hy(1))*(1/pi^2)* ...
%             (cos(pi*xpp(1:end-1,:))-cos(pi*xpp(2:end,:))).*...
%             (sin(pi*test(:,1:end-1))-sin(pi*test(:,2:end))) *t_factor2;
% 
% test = kron(ones(Nvy_in,1),x);
% test = reshape(test,Nx+1,Nvy_in);
% v_exact = 1/(hx(1)*hy(1))*(1/pi^2)* ...
%             (-sin(pi*test(1:end-1,:))+sin(pi*test(2:end,:))).*...
%             (cos(pi*ypp(:,1:end-1))-cos(pi*ypp(:,2:end)))*t_factor2;


% u_error = uh-u_exact(:);
% v_error = vh-v_exact(:);
% u_error_i(j,jj) = max(abs(u_error))
% v_error_i(j,jj) = max(abs(v_error))
% u_error_2(j,jj) = sqrt( sum(u_error.^2)/Nu )
% v_error_2(j,jj) = sqrt( sum(v_error.^2)/Nv )
% % 
% % k_exact = t_factor1;
% % k_error_i(j) = abs(k_exact(end)-k(end))
% % 
% 
p_exact      = 0.25*(cos(2*pi*xpp)+cos(2*pi*ypp))*t_factor1;
% p            = p - mean(p(:));
% p_error      = p(:)-p_exact(:);
% p_error_i(j,jj) = max(abs(p_error))
% p_error_2(j,jj) = sqrt( sum(p_error.^2)/Np)


%% calculate gamma, where e = C*(h^alfa) + K*(h^gamma)*(tau^beta)

% m = (jj-1)*3 + j
% if (m == 1)
%     usave = zeros(6,1);
%     vsave = zeros(6,1);
%     psave = zeros(6,1);
% end

if (jj==1)
    
    ufine = uh;
    vfine = vh;
    pfine = p-mean(p(:));
    
end
if (jj>1 || length(dt_list)==1)

    % temporal error obtained by subtracting spatial error (obtained with
    % ufine)
    u_error = uh-u_exact(:);
    v_error = vh-v_exact(:);
    ufine_error = ufine-u_exact(:);
    vfine_error = vfine-v_exact(:);
    
    p       = p - mean(p(:));
    p_error = p - p_exact(:);
    pfine_error = pfine-p_exact(:);
    
    
    % old way of calculating error:
    u_error_i(j,jj) = max(abs(uh-ufine))
    u_error_2(j,jj) = sqrt(sum( (uh-ufine).^2)/Nu)
    v_error_i(j,jj) = max(abs(vh-vfine))
    v_error_2(j,jj) = sqrt(sum( (vh-vfine).^2)/Nv)
    p_error_i(j,jj) = max(abs(p-pfine))
    p_error_2(j,jj) = sqrt(sum( (p-pfine).^2)/Np)
    
    % new way of calculating error, using exact solution
%     u_error_i(j,jj) = max(abs(u_error)) - max(abs(ufine_error))
%     v_error_i(j,jj) = max(abs(v_error)) - max(abs(vfine_error))
%     u_error_2(j,jj) = sqrt( sum(u_error.^2)/Nu ) - sqrt( sum(ufine_error.^2)/Nu)
%     v_error_2(j,jj) = sqrt( sum(v_error.^2)/Nv ) - sqrt( sum(vfine_error.^2)/Nv)
% 
%     p_error_i(j,jj) = max(abs(p_error)) - max(abs(pfine_error))
%     p_error_2(j,jj) = sqrt( sum(p_error.^2)/Np) - sqrt( sum(pfine_error.^2)/Np)
%     
    uspace_error(j) = max(abs(ufine_error));
    vspace_error(j) = max(abs(vfine_error));
    pspace_error(j) = max(abs(pfine_error));
end

if (jj==length(dt_list) && j==length(mesh_list))
    
    % spatial order
    disp('spatial order:')
    get_slope(1./mesh_list,abs(uspace_error))
    get_slope(1./mesh_list,abs(vspace_error))
    get_slope(1./mesh_list,abs(pspace_error))
    figure(2)
    loglog(1./mesh_list,abs(uspace_error),'bx-');
    hold on
    loglog(1./mesh_list,abs(vspace_error),'rx-');
    loglog(1./mesh_list,abs(pspace_error),'kx-');
    
    colors1 = {'kx-','rx-','bx-','gx-','mx-','yx-','cx-'};
    colors2 = {'ko--','ro--','bo--','go--','mo--','yo--','co--'};
    
    for k = 1:length(mesh_list)
        figure(1)
        loglog(dt_list(2:end),abs(u_error_i(k,2:end)),char(colors1(k)))
        hold on
%         figure(2)
        loglog(dt_list(2:end),abs(u_error_2(k,2:end)),char(colors2(k)))
%         hold on    
        
        % temporal order for each grid
        get_slope(dt_list(2:end),abs(u_error_i(k,2:end)),1,2)
        get_slope(dt_list(2:end),abs(u_error_2(k,2:end)),1,2)
    
    end
    
    legend('Linf N=5','L2 N=5','Linf N=10','L2 N=10','Linf N=20','L2 N=20',...
           'Linf N=40','L2 N=40');
    
    % average order reduction
    log(abs(u_error_i(2:end,3))./abs(u_error_i(1:end-1,3)))/log(2)
    log(abs(u_error_2(2:end,3))./abs(u_error_2(1:end-1,3)))/log(2)
    
end
% if (m == 6)
%     
%     dx1   = 2/mesh_list(1);
%     dx2   = 2/mesh_list(2);
%     
%     % spatial convergence for finest timestep, assuming negligible temporal
%     % error
%     alfa  = log ( usave(3) / usave(6) ) / log (dx1/dx2)
% 
%     % temporal convergence for both grids
%     beta  = log ( abs(usave(1)-usave(3)) / abs(usave(2)-usave(3)) ) / log(dt_list(1)/dt_list(2))
%     beta2 = log ( abs(usave(4)-usave(6)) / abs(usave(5)-usave(6)) ) / log(dt_list(1)/dt_list(2))            
% 
%     % this is K h_1^gamma
%     expr1 = abs(usave(1) - usave(3)) / (dt_list(1)^beta)
%     % this is K h_2^gamma
%     expr2 = abs(usave(4) - usave(6)) / (dt_list(1)^beta2)
%             
%     gamma = log( expr1/expr2) / log(dx1/dx2)   
%     
%     K     = expr1/(dx1^gamma)        
%     
% end