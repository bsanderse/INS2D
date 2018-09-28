%  Unsteady calculations:                                                   
%     - calculate convective and diffusive terms, and ustar and vstar (Ru
%     and Rv)
%     - calculate right-hand side vector of Poisson equation
%     - solve the Poisson equation for the pressure                                 
%     - update velocity field
%     - solve for k and epsilon

if (restart.load == 0)
    fprintf(fconv,'n            dt               t                res              maxdiv           umom             vmom             k\n');
    fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n',... 
                    n,dt,t,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
end


dp     = zeros(Np,1);
p      = zeros(Np,1);

dtn    = dt;
Newton = 1;

% zero block for total matrix
Z2     = spalloc(Np,Np,0);

% right hand side
f      = zeros(Nu+Nv+Np,1);

% while (abs(t)<=t_end-dt+eps)
while(n<=nt)
    
      n = n+1;    
           

      %% convective terms, u-component
      % c^n * u^(n+1), c=u
      uIux       = Iu_ux*uh+yIu_ux;                     % convective velocity, u_bar
      uAux       = Au_ux*uh+yAu_ux;
      C1         = Cux*spdiags(uIux,0,N1,N1);   
      C2         = Cux*spdiags(uAux,0,N1,N1)*Newton;
      Conv_ux_11 = C1*Au_ux + C2*Iu_ux;
      yConv_ux   = C1*uAux;
      
      % c^n * u^(n+1), c=v
      vIuy       = Iv_uy*vh+yIv_uy;                     % convective velocity, v_bar
      uAuy       = Au_uy*uh+yAu_uy;
      C1         = Cuy*spdiags(vIuy,0,N2,N2);   
      C2         = Cuy*spdiags(uAuy,0,N2,N2)*Newton;
      Conv_uy_11 = C1*Au_uy;
      Conv_uy_12 = C2*Iv_uy;
      yConv_uy   = C1*uAuy;
      

      %% convective terms, v-component
      % c^n * v^(n+1), c=u
      uIvx       = Iu_vx*uh+yIu_vx;                 % convective velocity, u_bar  
      vAvx       = Av_vx*vh+yAv_vx;
      C1         = Cvx*spdiags(uIvx,0,N3,N3);
      C2         = Cvx*spdiags(vAvx,0,N3,N3)*Newton;
      Conv_vx_21 = C2*Iu_vx;
      Conv_vx_22 = C1*Av_vx;
      yConv_vx   = C1*vAvx;
      
      % c^n * v^(n+1), c=v
      vIvy       = Iv_vy*vh+yIv_vy;                 % convective velocity, v_bar
      vAvy       = Av_vy*vh+yAv_vy;
      C1         = Cvy*spdiags(vIvy,0,N4,N4);
      C2         = Cvy*spdiags(vAvy,0,N4,N4)*Newton;       
      Conv_vy_22 = C1*Av_vy + C2*Iv_vy;
      yConv_vy   = C1*vAvy;
      
        
      %% diffusive terms

      % determine nu at different points
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
      
      
      Diffu_11   = Dux*( nu_ux * 2 * Su_ux) + Duy*( nu_uy * Su_uy);
      Diffu_12   = Duy*( nu_uy * Sv_uy);
      Diffv_21   = Dvx*( nu_vx * Su_vx);
      Diffv_22   = Dvx*( nu_vx * Sv_vx) + Dvy*( nu_vy * 2 * Sv_vy);

      yDiffu     = Dux*( nu_ux* 2 * (Su_ux*uh + ySu_ux)) + ...
                   Duy*( nu_uy* (Su_uy*uh + ySu_uy + Sv_uy*vh + ySv_uy));
      yDiffv     = Dvx*( nu_vx* (Su_vx*uh + ySu_vx + Sv_vx*vh + ySv_vx)) + ...
                   Dvy*( nu_vy* 2 * (Sv_vy*vh + ySv_vy));        

      CD_11      = spdiags(Omu,0,Nu,Nu)/dt + Conv_ux_11 + Conv_uy_11 - Diffu_11;
      CD_12      = Conv_uy_12 - Diffu_12;
      CD_21      = Conv_vx_21 - Diffv_21;
      CD_22      = spdiags(Omv,0,Nv,Nv)/dt  + Conv_vx_22 + Conv_vy_22 - Diffv_22;     
      CD         = [CD_11 CD_12; CD_21 CD_22];

      % right-hand side; this is the residual
      f(1:Nu)        = yConv_ux + yConv_uy - yDiffu - Fx + Gx*p + y_px;
      f(Nu+1:Nu+Nv)  = yConv_vx + yConv_vy - yDiffv - Fy + Gy*p + y_py;
          
         
      %% solve system for u,v,p
      % direct solver from Matlab
      f(Nu+Nv+1:end) = M*V + yM; 
      Z          = [CD G; M Z2];

      dq        = -Z\f;

      dV        = dq(1:Nu+Nv);
      dp        = dq(Nu+Nv+1:end);
      du        = dV(1:Nu);
      dv        = dV(Nu+1:end);     

      V         = V + dV;
      uh        = uh + du;
      vh        = vh + dv;
      p         = p + dp;
      
%       % pressure correction
%       dV        = -CD\f(1:Nu+Nv);
% 
%       du        = dV(1:Nu);
%       dv        = dV(Nu+1:end);     
% 
%       V         = V + dV;          
%       R         = V;
% 
%       pressure_correction;
% 
%       uh     = V(1:Nu);
%       vh     = V(Nu+1:end);
%       p      = p + dp;
              
    
      
      %% construct matrix for k-equation
      uIkx       = Iu_kx*uh+yIu_kx;
      vIky       = Iv_ky*vh+yIv_ky;
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
      
      nu_p       = Cmu*(kth.^2)./(eh + eps);
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
      
      % matrix
      Convk_k    = Ckx*uIkx_mat*Ak_kx + Cky*vIky_mat*Ak_ky;
      Diffk_k    = -Dkx*( nu_kx * Skx) - Dky*( nu_ky * Sky) + ...
                   -Dkx*( Cmu/sigmak * spdiags(2*kAkx./eAex .* kSkx,0,(Npx+1)*Npy,(Npx+1)*Npy))*Ak_kx + ...
                   -Dky*( Cmu/sigmak * spdiags(2*kAky./eAey .* kSky,0,Npx*(Npy+1),Npx*(Npy+1)))*Ak_ky;       
      Prodk_k    = -2*Cmu*spdiags( Omp.* kth./(eh+eps) .*prod,0,Np,Np);      
      Zk         = spdiags(Omp,0,Np,Np)/dt  + Convk_k + Diffk_k + ...
                   2*spdiags(Omp.*eh./(kth+eps),0,Np,Np);
% Zijlema: no Prodk_k, but dissipation rewritten:
%                   2*spdiags(Omp.*eh./(kth+eps),0,Np,Np);

      dkt        = -Zk\fk;      
      kth        = kth + dkt;
      
      %% construct matrix for e-equation
      
      eSex       = Sex*eh + ySex;
      eSey       = Sey*eh + ySey;
      nu_p       = Cmu*(kth.^2)./(eh + eps);
      % this uses the old kth:
      nu_ex      = spdiags(nu + nu_x/sigmae,0,(Npx+1)*Npy,(Npx+1)*Npy);
      nu_ey      = spdiags(nu + nu_y/sigmae,0,Npx*(Npy+1),Npx*(Npy+1));
      
      % residual
      yConve     = Ckx*(uIkx.*eAex) + Cky*(vIky.*eAey);
      yDiffe     = Dkx*( nu_ex * eSex ) + Dky*( nu_ey * eSey ); 
      yProde     = Ce1*(eh./(kth-dkt + eps)).*yProdk;      % Zijlema: use old kth, i.e. kth-dkt
      yDisse     = Ce2*(eh./(kth + eps)).*yDissk;      
      
      fe         = yConve - yProde + yDisse - yDiffe - Fe;
      
      % matrix
      Conve_e    = Ckx*uIkx_mat*Ae_ex + Cky*vIky_mat*Ae_ey;
      Disse_e    = 2*Ce2*spdiags(Omp.*eh./(kth+eps),0,Np,Np);
      Diffe_e    = -Dkx*( nu_ex * Sex) - Dky*( nu_ey * Sey) + ...
                    Dkx*( Cmu/sigmae * spdiags( (kAkx./eAex).^2 .* eSex,0,(Npx+1)*Npy,(Npx+1)*Npy))*Ae_ex + ...
                    Dky*( Cmu/sigmae * spdiags( (kAky./eAey).^2 .* eSey,0,Npx*(Npy+1),Npx*(Npy+1)))*Ae_ey;
      Ze         = spdiags(Omp,0,Np,Np)/dt + Conve_e + Disse_e  + Diffe_e;


      de         = -Ze\fe;      
      eh         = eh + de;      
      
      maxres(n)  = max([max(abs(f)) max(abs(fk)) max(abs(fe))]);

      % the velocities and pressure that are just computed are at 
      % the new time level t+dt:        
      t          = t + dt;
       
        
      process_iteration;

      % write convergence information to file  
      fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n',...                     
                    n,dt,t,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
                
      if (cw_output==1)
        fprintf(fcw,'%-10i %16.8e %16.8e %16.8e\n',n,dt,t,maxres(n));          
      end
                

end