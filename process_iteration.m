% process data from this iteration

% calculate mass, momentum and energy
[maxdiv(n), umom(n), vmom(n), k(n)] = check_conservation(V,t,options);
% note, for the ROM, we could also precompute:
% kinetic energy: 
% 0.5*V'*Om*V = 0.5*(B*R + Vbc)'*Om*(B*R + Vbc)
%    = 0.5*(R'*B'*Om*B*R + 2*Vbc'*Om*(B*R) + Vbc'*Om*Vbc)
%    = 0.5*(R'*R + (2*Vbc'*Om*B)*R + Vbc'*Om*Vbc)
% momentum:
% total mom = e*Om*V; (e=ones(NV,1))
%   = e*Om*(B*R+Vbc) = (e*Om*B)*R + e*Om*Vbc
% divergence:
% M*(B*R + Vbc) = (M*B)*R + M*Vbc


% residual (in Finite Volume form)
% for ke model residual also contains k and e terms and is computed
% in solver_unsteady_ke
if (options.rom.rom == 0)
    if (~strcmp(visc,'keps'))
        [maxres(n), ~, ~] = F(V,V,p,t,options,0);
    end
elseif (options.rom.rom == 1)
    % get ROM residual
    [maxres(n),~,~] = F_ROM(R,0,t,options,0);   
end

% statistics
%       if (statistics.write==1 && rem(n,statistics.n)==0)
%         write_statistics;
%       end

% postprocess energy spectrum along line
%       up    = reshape( Bup*(Au_ux * uh + yAu_ux), Npx, Npy);
%       vp    = reshape( Bvp*(Av_vy * vh + yAv_vy), Npx, Npy);
%       kp    = up.^2+vp.^2;
%       omega = Wv_vx*vh - Wu_uy*uh;
%       enstrophy = reshape(omega.^2,Npx,Npy);
%
%       if (rem(n,10)==0)
%
%         [Ek1(:,n/10), freq, ~]  = PSDfunc(kp(Npx/2,:),yp);
%         Ek2(:,n/10) = PSDfunc(kp(:,3*Npy/4),xp);
%         Ew1(:,n/10) = PSDfunc(enstrophy(Npx/2+1,:),yin);
%         Ew2(:,n/10) = PSDfunc(enstrophy(:,3*Npy/4),xin);
%       end


% change timestep based on operators
if (steady==0 && timestep.set==1 && rem(n,timestep.n) == 0)
    set_timestep;
end

% write restart file
if (restart.write==1 && rem(n,restart.n) == 0)
    write_restart;
end


%for real time plotting:
if (rtp.show==1 && rem(n-1,rtp.n) == 0)
    run(rtp.file);
    
    if (steady==0)
        fprintf(fcw,['t=' num2str(t) '\n']);
    end
    
    if (rtp.movie == 1)
        frame = getframe(gcf);
        writeVideo(writerObj,frame);
    end
    
    pause(0.01)
    
end

% store unsteady data in an array
if (steady==0 && save_unsteady == 1)
    uh_total(n,:) = V(1:options.grid.Nu);
    vh_total(n,:) = V(options.grid.Nu+1:end);
    p_total(n,:)  = p;
end

% write data to Tecplot file
if (tecplot.write==1 && rem(n,tecplot.n)==0)
    
    fprintf(fcw,'writing data to tecplot file... \n');
    Npx = options.grid.Npx;
    Npy = options.grid.Npy;

    [up,vp,qp] = get_velocity(V,t,options);
    
    pp  = reshape(p,Npx,Npy);
    %             Tp    = reshape(Tp,Nx,Ny);
    %             Tp    = zeros(Nx,Ny);
    Tp  = atand(vp./up); % local flow angle in degrees
    
    nzeros= filelen-length(num2str(n));
    n_new = num2str(10^nzeros);
    
    velocityTecplot2D([path_results '/' project n_new(2:end) num2str(n)], ...
        xp,yp,up,vp,pp,Tp);
    
end

% write convergence information to file
if (options.output.save_results == 1)
    if (steady == 0)
        fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n',...
            n,dt,t,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
    elseif (steady == 1)
        fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e \n',...
            n,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
    end
end