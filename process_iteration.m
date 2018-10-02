% process data from this iteration

      % calculate mass, momentum and energy
      [maxdiv(n), umom(n), vmom(n), k(n)] = check_conservation(uh,vh,t,options);
      
%       check_conservation;
            
      % residual (in Finite Volume form)
      % for ke model residual also contains k and e terms and is computed
      % in solver_unsteady_ke
      if (~strcmp(visc,'turbulent'))
          [maxres(n), ~, ~] = F(uh,vh,p,t,options);        
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
      if (rtp.show==1 && rem(n,rtp.n) == 0)
            switch rtp.type
                case {'vorticity'}
                    vorticity;
                case {'velocity'}
                    velocity;
                case {'quiver'}
                    velocity_vectors;
                case {'streamfunction'}
                    streamfunction;
                case {'pressure'}
                    pressure;
            end
             fprintf(fcw,['t=' num2str(t) '\n']);
             pause(0.1)
      end
        
      % write data to Tecplot file        
      if (tecplot.write==1 && rem(n,tecplot.n)==0)
          
            fprintf(fcw,'writing data to tecplot file... \n');
            up    = reshape( Bup*(Au_ux * uh + yAu_ux), Npx, Npy);
            vp    = reshape( Bvp*(Av_vy * vh + yAv_vy), Npx, Npy);
           
%             up = up - cos(t);
            
            pp    = reshape(p,Npx,Npy);
%             Tp    = reshape(Tp,Nx,Ny);
%             Tp    = zeros(Nx,Ny);
            Tp    = atand(vp./up); % local flow angle in degrees
            
            nzeros= filelen-length(num2str(n));
            n_new = num2str(10^nzeros);
           
            velocityTecplot2D([path_results '/' project n_new(2:end) num2str(n)], ...
                xp,yp,up,vp,pp,Tp);
           
      end