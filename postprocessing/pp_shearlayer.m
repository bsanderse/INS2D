%% compare velocity field with velocity field at fine mesh

% u_coarse = reshape(uh,Nux_in,Nuy_in);
% v_coarse = reshape(vh,Nvx_in,Nvy_in);
% xin_coarse = xin;
% yin_coarse = yin;
% xp_coarse  = xp;
% yp_coarse  = yp;
% Npx_coarse = Npx;

% load fine mesh data
% load_file = 'results/regularization/nomodel_Re500/N160_t8_nomodel.mat';
% load(load_file,'uh','vh','xin','yp','yin','xp','Npx','Npy');

% if (Npx~=Npx_coarse)
%     
%     % interpolate to coarse mesh positions
%     u_fc = interp2(yp',xin,reshape(uh,Npx,Npy),yp_coarse',xin_coarse);
%     v_fc = interp2(yin',xp,reshape(vh,Npx,Npy),yin_coarse',xp_coarse);
% 
%     error_u = max(abs(u_fc(:) - u_coarse(:)))
%     error_v = max(abs(v_fc(:) - v_coarse(:)))
% 
% else
%     error_u = max(abs(uh(:) - u_coarse(:)))
%     error_v = max(abs(vh(:) - v_coarse(:)))
% end

if (j==1)
    Nfine = Npx;
    ufine = uh;
    vfine = vh;
    xin_fine = xin;
    yin_fine = yin;
    xp_fine  = xp;
    yp_fine  = yp;    
    
else
    
    % interpolate to coarse mesh positions
    u_fc = interp2(yp_fine',xin_fine,reshape(ufine,Nfine,Nfine),yp',xin,'spline');
    v_fc = interp2(yin_fine',xp_fine,reshape(vfine,Nfine,Nfine),yin',xp,'spline');

    error_u(j) = max(abs(u_fc(:) - uh))
    error_v(j) = max(abs(v_fc(:) - vh))    
   
    
end

vorticity;
% pause