%  Unsteady calculations:                                                   
%     - calculate convective and diffusive terms, and ustar and vstar (Ru
%     and Rv)
%     - calculate right-hand side vector of Poisson equation
%     - solve the Poisson equation for the pressure                                 
%     - update velocity field

if (restart.load == 0)
    fprintf(fconv,'n            dt               t                res              maxdiv           umom             vmom             k\n');
    fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n',... 
                    n,dt,t,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));
end

if (method==2)
    cu     = uh;
    cv     = vh;
    convection;
    Cu_old = du2dx + duvdy;
    Cv_old = duvdx + dv2dy;
end

uh_old = uh;
vh_old = vh;
p_old  = p;
dp     = zeros(Np,1);

dtn    = dt;


% for methods that need a velocity field at n-1 the first time step
% use RK4 or FC2: one-leg, AB, CN, FC1, IM
% CN needs start-up for extrapolated Picard linearization
if (method==2 || method==4 || method==5 || ...
    method==71 || (method==61 && EP==1))
      n           = n+1;
      fprintf(fcw,['first time step method ' num2str(method_startup) '\n']);
      method_temp = method;
      method      = method_startup;

      % select time method
      select_time_method;
      
      method      = method_temp;
      
      % the velocities and pressure that are just computed are at 
      % the new time level t+dt:       
      t = t + dt;      
            
      % check residuals, conservation, write output files
      process_iteration;            
      
      % write convergence information to file  
      fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n',...                     
                    n,dt,t,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));   

end

eps = 1e-12;

% while (abs(t)<=(t_end-dt+eps))
while(n<=nt)

    % time reversal
%       if (n==nt/2+1)
%           fprintf(fcw,['reversing time at ' num2str(t) '\n']);
%           dt=-dt;
%       end

%%
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
%%
    
      n = n+1;    
           
      % perform one time step with the time integration method
      select_time_method;

      
      % the velocities and pressure that are just computed are at 
      % the new time level t+dt:       
      t = t + dt;
      
            
      % check residuals, conservation, write output files
      process_iteration;      
      

      % write convergence information to file  
      fprintf(fconv,'%-10i %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e \n',...                     
                    n,dt,t,maxres(n),maxdiv(n),umom(n),vmom(n),k(n));      
                
      
end