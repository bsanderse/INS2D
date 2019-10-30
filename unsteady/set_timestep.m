% estimate time step based on eigenvalues of operators, 
% using Gershgorin

% for explicit methods only
if (method==1 || method==2 || method==5 || method==81 || method==82)

    %% convective part
    Cu = Cux*spdiags(Iu_ux*uh+yIu_ux,0,N1,N1)*Au_ux + ...
         Cuy*spdiags(Iv_uy*vh+yIv_uy,0,N2,N2)*Au_uy;
    Cv = Cvx*spdiags(Iu_vx*uh+yIu_vx,0,N3,N3)*Av_vx + ...
         Cvy*spdiags(Iv_vy*vh+yIv_vy,0,N4,N4)*Av_vy; 

    test = spdiags(Omu_inv,0,Nu,Nu)*Cu;
    sum_conv_u = abs(test)*ones(Nu,1) - diag(abs(test)) - diag(test);
    test = spdiags(Omv_inv,0,Nv,Nv)*Cv;
    sum_conv_v = abs(test)*ones(Nv,1) - diag(abs(test)) - diag(test);
    labda_conv = max([max(sum_conv_u) max(sum_conv_v)]);

    %% diffusive part
    test = spdiags(Omu_inv,0,Nu,Nu)*Diffu;
    sum_diff_u = abs(test)*ones(Nu,1) - diag(abs(test)) - diag(test);
    test = spdiags(Omv_inv,0,Nv,Nv)*Diffv;
    sum_diff_v = abs(test)*ones(Nv,1) - diag(abs(test)) - diag(test);
    labda_diff = max([max(sum_diff_u) max(sum_diff_v)]);

    % based on max. value of stability region    
    if (method==5)
        labda_diff_max = 4*beta/(2*beta+1);
    elseif (method==1)
        labda_diff_max = 2;
    elseif (method==2)
        labda_diff_max = 1;
    elseif (method==82 || method==81)
        labda_diff_max = 2.78;
    end

    dt_diff = labda_diff_max/labda_diff;

    % based on max. value of stability region (not a very good indication
    % for the methods that do not include the imaginary axis)
    if (method==81 || method==82)
        labda_conv_max = 2*sqrt(2);
    elseif (method==11)
        labda_conv_max = sqrt(3);
    else
        labda_conv_max = 1;
    end
    dt_conv = labda_conv_max/labda_conv;

    dt = CFL*min(dt_conv,dt_diff);
    
end

% alternative: full eigenvalue analysis:

%         Cu = Cux*spdiags(Iu_ux*uh+yIu_ux,0,N1,N1)*Au_ux + ...
%              Cuy*spdiags(Iv_uy*vh+yIv_uy,0,N2,N2)*Au_uy;
%         Cv = Cvx*spdiags(Iu_vx*uh+yIu_vx,0,N3,N3)*Av_vx + ...
%              Cvy*spdiags(Iv_vy*vh+yIv_vy,0,N4,N4)*Av_vy;
%
%         Eu = eig(full(spdiags(Omu_inv,0,Nu,Nu)*Cu));
%         Ev = eig(full(spdiags(Omv_inv,0,Nv,Nv)*Cv));
%         max_eig(n) = max(abs([Eu;Ev]))*dt;
%         dtn = dt;
%         set_timestep;
%         dt  = dtn;
%         max_eig_G(n) = labda_conv*dt;