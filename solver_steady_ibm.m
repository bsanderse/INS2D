% solve the entire saddlepoint system arising from the steady Navier-Stokes
% equations with linearization of the convective terms


% right hand side
f = zeros(Nu+Nv+Np,1);
p = zeros(Np,1);

Newton = 0;

% remove rows of M corresponding to inside and interface pressure points
% sorting is not strictly necessary, but increases matrix readability and
% multiplication

% divergence matrix for interior and interface points
M_rhs              = M(sort(index_p),:);   
% divergence matrix for outer points
M_outer            = M;
M_outer(index_p,:) = [];   
yM_outer           = yM;
yM_outer(index_p)  = [];

% change size right-hand side vector
f(Nu+Nv+index_p)   = [];

% pressure in outer points
p_outer            = p(:);     
p_outer(index_p)   = [];


% gradient operator to be added to right hand side
G_rhs              = -M_rhs';
% new gradient operator based on adapted divergence operator
G_outer            = -M_outer';
Gx_outer           = G_outer(1:Nu,:);
Gy_outer           = G_outer(Nu+1:Nu+Nv,:);

% size of zero block in saddle-point matrix
Z2                 = spalloc(size(M_outer,1),size(M_outer,1),0);

maxres(1)=10*accuracy;
while ( maxres(n) > accuracy)

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
      vAvx       = Av_vx*vh + yAv_vx;
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
      
      
      %% construct matrix (saddlepoint structure)
      CD_11      = - Diffu + Conv_ux_11 + Conv_uy_11;
      CD_12      = Conv_uy_12;
      CD_21      = Conv_vx_21;
      CD_22      = - Diffv + Conv_vx_22 + Conv_vy_22;     
      CD         = [CD_11 CD_12; CD_21 CD_22];
           
      % right-hand side; this is -1*residual
      f(1:Nu)        = Diffu*uh + yDiffu - yConv_ux - yConv_uy + ...
                       Fx - Gx_outer*p_outer - y_px;
      f(Nu+1:Nu+Nv)  = Diffv*vh + yDiffv - yConv_vx - yConv_vy + ...
                       Fy - Gy_outer*p_outer - y_py;
      f(Nu+Nv+1:end) = -M_outer*V - yM_outer;
      

      %% adapt rows that involve interface points
      % matrix sizes unaltered  
      
      % pressure
      rhs_p    = set_interpolation_rhs(xbi_p,ybi_p,xip_p,yip_p,xp,yp,interface_p,reshape(p,Npx,Npy),pB);         
      for qq=1:n_interface_p

          pGC(qq)           = beta_p(qq,:)*rhs_p(qq,:)';
       
          % set off-diagonals to interpolation value
%           [dum1 dum2 val]   = find(bilin_points_p(qq,:));
%           G_outer(ind,val)  = beta_p(qq,dum2);
       
      end

      
      % add pGC contribution to rhs
      p_inner(interface_p_inner) = pGC(interface_p_order);      
      f(1:Nu+Nv)                 = f(1:Nu+Nv) - G_rhs*p_inner;
      
      
      % u-velocity
      rhs_u                 = set_interpolation_rhs(xbi_u,ybi_u,...
                               xip_u,yip_u,xin,yp,interface_u,u,uB);
      CD(interface_u_ind,:) = 0;          
      for qq=1:n_interface_u


          ind               = interface_u_ind(qq);
         
          % set off-diagonals to interpolation value
          [dum1 dum2 val]   = find(bilin_points_u(qq,:));
          CD(ind,val)       = beta_u(qq,dum2);
          
          % adapt rhs
          f(ind)            = 2*uB - uh(ind) - beta_u(qq,:)*rhs_u(qq,:)';
          
      end

      
      %v-velocity
      rhs_v                 = set_interpolation_rhs(xbi_v,ybi_v,...
                               xip_v,yip_v,xp,yin,interface_v,v,vB);         
      CD(interface_v_ind+Nu,:) = 0;      
      for qq=1:n_interface_v

          ind               = interface_v_ind(qq);
         
          % set off-diagonals to interpolation value
          [dum1 dum2 val]   = find(bilin_points_v(qq,:));
          CD(ind+Nu,val+Nu) = beta_v(qq,dum2);
          
          % adapt rhs
          f(ind+Nu)         = 2*vB - vh(ind) - beta_v(qq,:)*rhs_v(qq,:)';          

      end

      % set diagonal of interface points to 1     
      CD = CD + CD_interface;
      
      
      %% adapt rows that involve inside points (not interface)
      % set row to zero
      CD(index_u,:)    = 0;
      f(index_u)       = uB - uh(index_u);
          
      CD(index_v+Nu,:) = 0;     
      f(index_v+Nu)    = vB - vh(index_v);      
      
      % set diagonal of inside points to 1
      CD               = CD + CD_inside;
      
      
      %% total matrix
      Z         = [CD G_outer; M_outer Z2];

      
      %% solve with direct solver from Matlab
      dq        = Z\f;
     
      dV        = dq(1:Nu+Nv);
      dp        = dq(Nu+Nv+1:end);
      du        = dV(1:Nu);
      dv        = dV(Nu+1:end);     
            
      V         = V + dV;
      uh        = uh + du;
      vh        = vh + dv;
      p_outer   = p_outer + dp;      
      p(index_p_outer) = p_outer;
      p(sort(index_p)) = p_inner;
      
      u         = reshape(uh,Nux_in,Nuy_in);
      v         = reshape(vh,Nvx_in,Nvy_in);

      % check residuals, conservation, write output files    
      process_iteration;
      maxres(n) = max(abs([resu.*(1-inside_incl_interface_u(:));...
                           resv.*(1-inside_incl_interface_v(:))]));
      
      maxres(n)

      % write convergence information to file  
      fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e \n',...
                    n,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
      
 
end