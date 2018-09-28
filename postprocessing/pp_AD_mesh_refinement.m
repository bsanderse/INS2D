% process velocity profiles 

if (j==1) % finest result
    
    ux_fine = ux; % @ xin
    uy_fine = uy; % @ yp
    vy_fine = vy; % @ yin
    px_fine = px; % @ xp
    py_fine = py; % @ yp
    yin_fine = yin;
    yp_fine = yp; 
    xin_fine = xin;
    xp_fine = xp;
    
    k_fine  = k(end);
    
else
    % calculate errors
    
    k_error(j-1) = abs(k_fine-k(end))
    
    % L2 errors
    ux_error_2(j-1) = sqrt((1/length(xin))*sum((interp1(xin_fine,ux_fine,xin,'linear') - ux).^2))
    uy_error_2(j-1) = sqrt((1/length(yp))*sum((interp1(yp_fine,uy_fine',yp,'linear') - uy').^2));
    px_error_2(j-1) = sqrt((1/length(xp))*sum((interp1(xp_fine,px_fine,xp,'linear') - px).^2));
    py_error_2(j-1) = sqrt((1/length(yp))*sum((interp1(yp_fine,py_fine',yp,'linear') - py').^2));
    vy_error_2(j-1) = sqrt((1/length(yin))*sum((interp1(yin_fine,vy_fine',yin,'linear') - vy').^2));

    
    % L_inf errors
    ux_error_i(j-1) = max(abs(interp1(xin_fine,ux_fine,xin,'linear') - ux));
    uy_error_i(j-1) = max(abs(interp1(yp_fine,uy_fine',yp,'linear') - uy'));
    px_error_i(j-1) = max(abs(interp1(xp_fine,px_fine,xp,'linear') - px));
    py_error_i(j-1) = max(abs(interp1(yp_fine,py_fine',yp,'linear') - py'));
    vy_error_i(j-1) = max(abs(interp1(yin_fine,vy_fine',yin,'linear') - vy'));
    
end
  
%%
if (j==length(mesh_list))
    
    figure
    loglog(1./mesh_list(2:end),k_error,'mx-')
    hold on
    loglog(1./mesh_list(2:end),ux_error_2,'bx-')
    loglog(1./mesh_list(2:end),uy_error_2,'bx--')
    loglog(1./mesh_list(2:end),vy_error_2,'rx-')
    loglog(1./mesh_list(2:end),px_error_2,'kx-')
    loglog(1./mesh_list(2:end),py_error_2,'kx--')
    
    loglog(1./mesh_list(2:end),ux_error_i,'bo-')
    loglog(1./mesh_list(2:end),uy_error_i,'bo--')
    loglog(1./mesh_list(2:end),vy_error_i,'ro-')
    loglog(1./mesh_list(2:end),px_error_i,'ko-')
    loglog(1./mesh_list(2:end),py_error_i,'ko--')
end