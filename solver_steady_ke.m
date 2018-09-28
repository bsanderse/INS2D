% solve the entire saddlepoint system arising from the steady Navier-Stokes
% equations with linearization of the convective terms

% zero block for total matrix
Z2 = spalloc(Np,Np,0);

% right hand side
f = ones(Nu+Nv+3*Np,1);

Newton = 0;
eps    = 1e-14; % prevent division by zero

      
while ( max(abs(f)) > accuracy)

      % switch to Newton after nPicard steps
      if (strcmp(linearization,'Newton') && n>nPicard)
          Newton = 1;
      end

      n = n+1; 
      
      % to check if the convection matrix C is skew-symmetric, look at
%       Cu = Cux*spdiags(Iu_ux*uh+yIu_ux,0,N1,N1)*Au_ux + Cuy*spdiags(Iv_uy*vh+yIv_uy,0,N2,N2)*Au_uy;
%       Cv = Cvx*spdiags(Iu_vx*uh+yIu_vx,0,N3,N3)*Av_vx + Cvy*spdiags(Iv_vy*vh+yIv_vy,0,N4,N4)*Av_vy;
%       max2d(abs(Cu+Cu'))
%       max2d(abs(Cv+Cv'))
     
      
      %% convective terms, u-component
      % c^n * u^(n+1), c=u
      uIux       = Iu_ux*uh + yIu_ux;                     % convective velocity, u_bar
      uAux       = Au_ux*uh + yAu_ux;
      C1         = Cux*spdiags(uIux,0,N1,N1);   
      C2         = Cux*spdiags(uAux,0,N1,N1)*Newton;
      Conv_ux_11 = C1*Au_ux + C2*Iu_ux;
      yConv_ux   = C1*uAux;
      
      % c^n * u^(n+1), c=v
      vIuy       = Iv_uy*vh + yIv_uy;                     % convective velocity, v_bar
      uAuy       = Au_uy*uh + yAu_uy;
      C1         = Cuy*spdiags(vIuy,0,N2,N2);   
      C2         = Cuy*spdiags(uAuy,0,N2,N2)*Newton;
      Conv_uy_11 = C1*Au_uy;
      Conv_uy_12 = C2*Iv_uy;
      yConv_uy   = C1*uAuy;
      

      %% convective terms, v-component
      % c^n * v^(n+1), c=u
      uIvx       = Iu_vx*uh + yIu_vx;                 % convective velocity, u_bar  
      vAvx       = Av_vx*vh + yAv_vx;
      C1         = Cvx*spdiags(uIvx,0,N3,N3);
      C2         = Cvx*spdiags(vAvx,0,N3,N3)*Newton;
      Conv_vx_21 = C2*Iu_vx;
      Conv_vx_22 = C1*Av_vx;
      yConv_vx   = C1*vAvx;
      
      % c^n * v^(n+1), c=v
      vIvy       = Iv_vy*vh + yIv_vy;                 % convective velocity, v_bar
      vAvy       = Av_vy*vh + yAv_vy;
      C1         = Cvy*spdiags(vIvy,0,N4,N4);
      C2         = Cvy*spdiags(vAvy,0,N4,N4)*Newton;       
      Conv_vy_22 = C1*Av_vy + C2*Iv_vy;
      yConv_vy   = C1*vAvy;
      
      
      %% diffusive terms
      
      % average k and e to ux, uy, vx, vy positions
      kAux       = Ak_ux*kth + yAk_ux;
      kAuy       = Ak_uy*kth + yAk_uy;
      kAvx       = Ak_vx*kth + yAk_vx;
      kAvy       = Ak_vy*kth + yAk_vy;
      eAux       = Ae_ux*eh  + yAe_ux + eps;
      eAuy       = Ae_uy*eh  + yAe_uy + eps;
      eAvx       = Ae_vx*eh  + yAe_vx + eps;
      eAvy       = Ae_vy*eh  + yAe_vy + eps;      
      nu_ux      = spdiags( nu + Cmu* kAux.^2 ./ eAux,0,N1,N1);
      nu_uy      = spdiags( nu + Cmu* kAuy.^2 ./ eAuy,0,N2,N2);
      nu_vx      = spdiags( nu + Cmu* kAvx.^2 ./ eAvx,0,N3,N3);
      nu_vy      = spdiags( nu + Cmu* kAvy.^2 ./ eAvy,0,N4,N4);
                  
      % difference u and v to ux, uy, vx, vy positions
      uSux       = Su_ux*uh + ySu_ux;
      uSuy       = Su_uy*uh + ySu_uy + Sv_uy*vh + ySv_uy;
      vSvx       = Su_vx*uh + ySu_vx + Sv_vx*vh + ySv_vx;
      vSvy       = Sv_vy*vh + ySv_vy;
                     
      % residual
      yDiffu     = Dux*( nu_ux* 2 * uSux ) + Duy*( nu_uy* uSuy );
      yDiffv     = Dvx*( nu_vx* vSvx ) + Dvy*( nu_vy* 2 * vSvy);
      fu         = yConv_ux + yConv_uy - yDiffu - Fx + Gx*p + y_px;
      fv         = yConv_vx + yConv_vy - yDiffv - Fy + Gy*p + y_py;       
      
      % matrices
      Diffu_11   = Dux*( nu_ux * 2 * Su_ux) + Duy*( nu_uy * Su_uy);
      Diffu_12   = Duy*( nu_uy * Sv_uy);
      Diffv_21   = Dvx*( nu_vx * Su_vx);
      Diffv_22   = Dvx*( nu_vx * Sv_vx) + Dvy*( nu_vy * 2 * Sv_vy);
      
      Zuu        = Conv_ux_11 + Conv_uy_11 - Diffu_11;
      Zuv        = Conv_uy_12 - Diffu_12;
      Zvu        = Conv_vx_21 - Diffv_21;
      Zvv        = Conv_vx_22 + Conv_vy_22 - Diffv_22;   
      
      diag1      = 2*Cmu* (kAux./eAux) * 2 .* uSux;
      diag2      = 2*Cmu* (kAuy./eAuy) .* uSuy;
      Zuk        = - Dux*spdiags(diag1,0,N1,N1)*Ak_ux + ...
                   - Duy*spdiags(diag2,0,N2,N2)*Ak_uy;
      diag1      = Cmu* (kAux ./ eAux).^2 * 2 .* uSux;
      diag2      = Cmu* (kAuy ./ eAuy).^2 .* uSuy;
      Zue        = Dux*spdiags(diag1,0,N1,N1)*Ae_ux + ...
                   Duy*spdiags(diag2,0,N2,N2)*Ae_uy;
      
      diag1      = 2*Cmu* (kAvx./eAvx) .* vSvx;
      diag2      = 2*Cmu* (kAvy./eAvy) * 2 .* vSvy;
      Zvk        = - Dvx*spdiags(diag1,0,N3,N3)*Ak_vx + ...
                   - Dvy*spdiags(diag2,0,N4,N4)*Ak_vy;
      diag1      = Cmu* (kAvx ./ eAvx).^2 .* vSvx;
      diag2      = Cmu* (kAvy ./ eAvy).^2 * 2 .* vSvy;
      Zve        = Dvx*spdiags(diag1,0,N3,N3)*Ae_vx + ...
                   Dvy*spdiags(diag2,0,N4,N4)*Ae_vy;                 

      
      %% construct matrix for k-equation
      uIkx       = Iu_kx*uh + yIu_kx;
      vIky       = Iv_ky*vh + yIv_ky;
      uIkx_mat   = spdiags(uIkx,0,(Npx+1)*Npy,(Npx+1)*Npy);
      vIky_mat   = spdiags(vIky,0,Npx*(Npy+1),Npx*(Npy+1));
      
      kAkx       = Ak_kx*kth + yAk_kx;
      kAky       = Ak_ky*kth + yAk_ky;
      eAex       = Ae_ex*eh  + yAe_ex + eps;
      eAey       = Ae_ey*eh  + yAe_ey + eps;
      
      uCkx       = Cux_k*uh + yCux_k;
      uCky       = Cuy_k*(Auy_k*uh+yAuy_k) + yCuy_k;
      vCkx       = Cvx_k*(Avx_k*vh+yAvx_k) + yCvx_k;
      vCky       = Cvy_k*vh + yCvy_k;
      
      nu_p       = Cmu * kth.^2  ./(eh + eps);
      nu_x       = Cmu * kAkx.^2 ./ eAex; 
      nu_y       = Cmu * kAky.^2 ./ eAey;
      nu_kx      = spdiags(nu + nu_x/sigmak,0,(Npx+1)*Npy,(Npx+1)*Npy);
      nu_ky      = spdiags(nu + nu_y/sigmak,0,Npx*(Npy+1),Npx*(Npy+1));

      prod       = 2*uCkx.^2 + 2*vCky.^2 + (uCky + vCkx).^2;

      kSkx       = Skx*kth + ySkx;
      kSky       = Sky*kth + ySky;
      
      % residual
      yConvk     = Ckx*(uIkx.*kAkx) + Cky*(vIky.*kAky);      
      yDiffk     = Dkx*( nu_kx * kSkx ) + Dky*( nu_ky * kSky ); 
    
      yProdk     = Omp.*nu_p.*prod;
      yDissk     = Omp.*eh;
      
      fk         = yConvk - yProdk + yDissk - yDiffk - Fk;
      
      % matrices
      Convk_u    = Ckx*spdiags(kAkx,0,(Npx+1)*Npy,(Npx+1)*Npy)*Iu_kx;
      Prodk_u    = -4*spdiags(Omp.* nu_p .* uCkx,0,Np,Np) * Cux_k + ...
                   -2*spdiags(Omp.* nu_p .* (uCky + vCkx),0,Np,Np) * Cuy_k*Auy_k;               
      Zku        = Convk_u + Prodk_u;
                   
      Convk_v    = Cky*spdiags(kAky,0,Npx*(Npy+1),Npx*(Npy+1))*Iv_ky;
      Prodk_v    = -2*spdiags(Omp.* nu_p .* (uCky + vCkx),0,Np,Np) * Cvx_k*Avx_k + ...
                   -4*spdiags(Omp.* nu_p .* vCky,0,Np,Np) * Cvy_k;
      Zkv        = Convk_v + Prodk_v;
               
      Convk_k    = Ckx*uIkx_mat*Ak_kx + Cky*vIky_mat*Ak_ky;
      Diffk_k    = -Dkx*( nu_kx * Skx) - Dky*( nu_ky * Sky) + ...
                   -Dkx*( Cmu/sigmak * spdiags(2*kAkx./eAex .* kSkx,0,(Npx+1)*Npy,(Npx+1)*Npy))*Ak_kx + ...
                   -Dky*( Cmu/sigmak * spdiags(2*kAky./eAey .* kSky,0,Npx*(Npy+1),Npx*(Npy+1)))*Ak_ky;       
      Prodk_k    = -2*Cmu*spdiags( Omp.* kth./(eh+eps) .*prod,0,Np,Np);
      Zkk        = Convk_k + Diffk_k + Prodk_k;

% Zijlema: no Prodk_k, but dissipation rewritten: (this removes the quadratic
% convergence of Newton!)
%                   2*spdiags(Omp.*eh./(kth+eps),0,Np,Np);
      
      Prodk_e    = spdiags(Omp.* (nu_p./(eh+eps).*prod + 1),0,Np,Np);
      Diffk_e    = Dkx*( Cmu/sigmak * spdiags( (kAkx./eAex).^2 .* kSkx,0,(Npx+1)*Npy,(Npx+1)*Npy))*Ae_ex + ...
                   Dky*( Cmu/sigmak * spdiags( (kAky./eAey).^2 .* kSky,0,Npx*(Npy+1),Npx*(Npy+1)))*Ae_ey;       
      Zke        = Prodk_e + Diffk_e;
                   
          
      
      %% construct matrix for e-equation
           
      eSex       = Sex*eh + ySex;
      eSey       = Sey*eh + ySey;
      nu_ex      = spdiags(nu + nu_x/sigmae,0,(Npx+1)*Npy,(Npx+1)*Npy);
      nu_ey      = spdiags(nu + nu_y/sigmae,0,Npx*(Npy+1),Npx*(Npy+1));
      
      % residual
      yConve     = Ckx*(uIkx.*eAex) + Cky*(vIky.*eAey);
      yDiffe     = Dkx*( nu_ex * eSex ) + Dky*( nu_ey * eSey ); 
      
      yProde     = Ce1*(eh./(kth + eps)).*yProdk;    
      yDisse     = Ce2*(eh./(kth + eps)).*yDissk;
            
      fe         = yConve - yProde + yDisse - yDiffe - Fe;
      
      % matrices    
      Conve_u    = Ckx*spdiags(eAex,0,(Npx+1)*Npy,(Npx+1)*Npy)*Iu_kx;
      Prode_u    = -4*Ce1*Cmu* spdiags(Omp.*uCkx.*kth,0,Np,Np) * Cux_k + ...
                   -2*Ce1*Cmu* spdiags(Omp.*(uCky + vCkx).*kth,0,Np,Np) * Cuy_k*Auy_k;
      Zeu        = Conve_u + Prode_u;
      
      Conve_v    = Cky*spdiags(eAey,0,Npx*(Npy+1),Npx*(Npy+1))*Iv_ky;
      Prode_v    = -2*Ce1*Cmu* spdiags(Omp.*(uCky + vCkx).*kth,0,Np,Np) * Cvx_k*Avx_k + ...
                   -4*Ce1*Cmu* spdiags(Omp.*vCky.*kth,0,Np,Np) * Cvy_k;
      Zev        = Conve_v + Prode_v;

      Prode_k    = -Ce1*Cmu*spdiags(Omp.*prod,0,Np,Np);
      Disse_k    = -Ce2*spdiags(Omp.*(eh./(kth + eps)).^2,0,Np,Np);
      Diffe_k    = -Dkx*( Cmu/sigmae * spdiags(2*kAkx./eAex .* eSex,0,(Npx+1)*Npy,(Npx+1)*Npy))*Ae_ex + ...
                   -Dky*( Cmu/sigmae * spdiags(2*kAky./eAey .* eSey,0,Npx*(Npy+1),Npx*(Npy+1)))*Ae_ey;
      Zek        = Prode_k + Disse_k + Diffe_k;    
                  
      Conve_e    = Ckx*uIkx_mat*Ae_ex + Cky*vIky_mat*Ae_ey;
      Disse_e    = 2*Ce2*spdiags(Omp.*eh./(kth+eps),0,Np,Np);
      Diffe_e    = -Dkx*( nu_ex * Sex) - Dky*( nu_ey * Sey) + ...
                    Dkx*( Cmu/sigmae * spdiags( (kAkx./eAex).^2 .* eSex,0,(Npx+1)*Npy,(Npx+1)*Npy))*Ae_ex + ...
                    Dky*( Cmu/sigmae * spdiags( (kAky./eAey).^2 .* eSey,0,Npx*(Npy+1),Npx*(Npy+1)))*Ae_ey;
      Zee        = Conve_e + Disse_e  + Diffe_e;
       


      %% construct matrix (saddlepoint structure)   
           
      Z          = [Zuu Zuv Gx Zuk Zue; ...
                    Zvu Zvv Gy Zvk Zve; ...
                    M       Z2 Z2  Z2;  ...
                    Zku Zkv Z2 Zkk Zke; ...
                    Zeu Zev Z2 Zek Zee];

      % right-hand side; this is the residual
      f(1:Nu)                  = fu;
      f(Nu+1:Nu+Nv)            = fv;
      f(Nu+Nv+1:Nu+Nv+Np)      = M*V + yM;
      f(Nu+Nv+Np+1:Nu+Nv+2*Np) = fk;
      f(Nu+Nv+2*Np+1:end)      = fe;
      

      %% solve with direct solver from Matlab      
      dq         = -Z\f;
     
      dV         = dq(1:Nu+Nv);
      dp         = dq(Nu+Nv+1:Nu+Nv+Np);
      dkt        = dq(Nu+Nv+Np+1:Nu+Nv+2*Np);
      de         = dq(Nu+Nv+2*Np+1:end);
      du         = dV(1:Nu);
      dv         = dV(Nu+1:end);     
            
      V          = V + dV;
      uh         = uh + du;
      vh         = vh + dv;
      p          = p + dp;
      eh         = eh + de;
      kth        = kth + dkt;

          
      % calculate mass, momentum and energy      
      check_conservation;    
      
      % right-hand side is directly the residual
      max_res(n) = max(abs(f));
      
      
      %% output during running
      if (n/numchecks == floor(n/numchecks)) 

        fprintf('iteration %3i, residual %0.8g\n',n,max_res(n));

        % write convergence information to temporary file
        tab = sprintf('\t');  
        fwrite(ftmp,[num2str(n) tab num2str(max_res(n)) nl]);

        % track CPU time
        time(n) = toc;        
        
        
        %for real time plotting:
        if rtp==1
            switch rtptype
                case {'vorticity'}
                    vorticity;
                case {'quiver'}
                    velocity_vectors;
                case {'velocity'}
                    velocity;
                case {'streamfunction'}
                    streamfunction;
                case {'pressure'}
                    pressure;
                case {'ke'}
                    ke;
            end
        end
        
        % write data to Tecplot file        
        if movie==1
            Bmap  = kron(speye(Nuy_in),BMx);
            up    = reshape( Bmap*(Au_ux * uh + yAu_ux), Npx, Npy);

            Bmap  = kron(BMy,speye(Nvx_in));
            vp    = reshape( Bmap*(Av_vy * vh + yAv_vy), Npx, Npy);            

            pp    = reshape(p,Nx,Ny);
            Tp    = zeros(Nx,Ny);
            
            nzeros= filelen-length(num2str(n));
            n_new = num2str(10^nzeros);
           
            velocityTecplot2D([tecpath '/' project n_new(2:end) num2str(n)], ...
                xp,yp,up,vp,pp,Tp);
           
        end
        
      end

end