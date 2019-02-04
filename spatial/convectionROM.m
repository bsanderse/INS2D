function [convu, convv, Jacu, Jacv] = convectionROM(ru,rv,t,options,getJacobian)
% evaluate convective terms and, optionally, Jacobians

order4     = options.discretization.order4;
regularize = options.case.regularize;

M = options.rom.M;
Bu = options.rom.Bu;
Bv = options.rom.Bv;

Cux = options.discretization.Cux;
Cuy = options.discretization.Cuy;
Cvx = options.discretization.Cvx;
Cvy = options.discretization.Cvy;

Au_ux = options.discretization.Au_ux;
Au_uy = options.discretization.Au_uy;
Av_vx = options.discretization.Av_vx;
Av_vy = options.discretization.Av_vy;

yAu_ux = options.discretization.yAu_ux;
yAu_uy = options.discretization.yAu_uy;
yAv_vx = options.discretization.yAv_vx;
yAv_vy = options.discretization.yAv_vy;

Iu_ux = options.discretization.Iu_ux;
Iv_uy = options.discretization.Iv_uy;
Iu_vx = options.discretization.Iu_vx;
Iv_vy = options.discretization.Iv_vy;

yIu_ux = options.discretization.yIu_ux;
yIv_uy = options.discretization.yIv_uy;
yIu_vx = options.discretization.yIu_vx;
yIv_vy = options.discretization.yIv_vy;

Nu = options.grid.Nu;
Nv = options.grid.Nv;

Jacu = spalloc(Nu,Nu+Nv,0);
Jacv = spalloc(Nv,Nu+Nv,0);





%%

% assumes uh, vh, cu, cv
% normally cu = uh, cv = vh

if (order4==0)
    
    %% no regularization
    if (regularize == 0)
        
        % precompute matrices
%         if (
        N1 = size(Iu_ux,1);
        N2 = size(Iv_uy,1);
        N3 = size(Iu_vx,1);
        N4 = size(Iv_vy,1);        
        Chat_ux = zeros(M,M,M); % this poses severe memory requirements (M^3)!       
        yhat_ux = zeros(M,M); 
        Chat_uy = zeros(M,M,M);
        yhat_uy = zeros(M,M);         
        Chat_vx = zeros(M,M,M); % this poses severe memory requirements (M^3)!       
        yhat_vx = zeros(M,M); 
        Chat_vy = zeros(M,M,M);
        yhat_vy = zeros(M,M);     
        
        % note Bu = Nu*m

        % preprocessing:
        % Bu'*du2dx  
        % = Bu'*Cux*diag(Iu_ux*uh+yIu_ux)*(Au_ux*uh+yAu_ux);
        % = Bu'*Cux*diag(Iu_ux*Bu*ru+yIu_ux)*(Au_ux*Bu*ru+yAu_ux)
        % = sum_i Bu'*Cux*diag(Iu_ux*Bu_i*ru_i+yIu_ux)*Au_ux*Bu*ru + 
        %         Bu'*Cux*diag(Iu_ux*Bu_i*ru_i+yIu_ux)*yAu_ux
        % =   (sum_i [ru_i (Bu'*Cux*diag(Iu_ux*Bu_i)] + Bu'*Cux*diag(yIu_ux))*Au_ux*Bu*ru
        %   + (sum_i [ru_i (Bu'*Cux*diag(yIu_ux)] + Bu'*Cux*diag(yIu_ux))*yAu_ux
        % = sum_i [ru_i (Bu'*Cux*diag(Iu_ux*Bu_i)]*Au_ux*Bu*ru + 
        %   sum_i [ru_i (Bu'*Cux*diag(Iu_ux*Bu_i)])*yAu_ux + 
        %   Bu'*Cux*diag(yIu_ux)*Au_ux*Bu*ru + 
        %   Bu'*Cux*diag(yIu_ux)*yAu_ux

        % similarly
        % Bu'*duvdy  
        % = Bu'*Cuy*diag(Iv_uy*vh+yIv_uy)*(Au_uy*uh+yAu_uy);
        % = Bu'*Cuy*diag(Iu_ux*Bv*rv+yIv_uy)*(Au_uy*Bu*ru+yAu_uy)
        % = sum_i Bu'*Cuy*diag(Iu_ux*Bv_i*rv_i+yIv_uy)*Au_uy*Bu*ru + 
        %         Bu'*Cuy*diag(Iu_ux*Bv_i*rv_i+yIv_uy)*yAu_uy
        % =   (sum_i [rv_i (Bu'*Cuy*diag(Iv_uy*Bv_i)] + Bu'*Cuy*diag(yIv_uy))*Au_uy*Bu*ru
        %   + (sum_i [rv_i (Bu'*Cuy*diag(yIv_uy)] + Bu'*Cuy*diag(yIv_uy))*yAu_uy
        % = sum_i [rv_i (Bu'*Cuy*diag(Iv_uy*Bv_i)]*Au_uy*Bu*ru + 
        %   sum_i [rv_i (Bu'*Cuy*diag(Iv_uy*Bv_i)])*yAu_uy + 
        %   Bu'*Cuy*diag(yIv_uy)*Au_uy*Bu*ru + 
        %   Bu'*Cuy*diag(yIv_uy)*yAu_uy        
        
        B1 = Bu'*Cux;
        B2 = Bu'*Cuy;
        B3 = Bv'*Cvx;
        B4 = Bv'*Cvy;        
        for i=1:M
            
            Chat_ux(i,:,:) = B1*spdiags(Iu_ux*Bu(:,i),0,N1,N1)*Au_ux*Bu;
            yhat_ux(i,:)   = B1*spdiags(Iu_ux*Bu(:,i),0,N1,N1)*yAu_ux;
            Chat_uy(i,:,:) = B2*spdiags(Iv_uy*Bv(:,i),0,N2,N2)*Au_uy*Bu;
            yhat_uy(i,:)   = B2*spdiags(Iv_uy*Bv(:,i),0,N2,N2)*yAu_uy;
            
            Chat_vx(i,:,:) = B3*spdiags(Iu_vx*Bu(:,i),0,N3,N3)*Av_vx*Bv;
            yhat_vx(i,:)   = B3*spdiags(Iu_vx*Bu(:,i),0,N3,N3)*yAv_vx;
            Chat_vy(i,:,:) = B4*spdiags(Iv_vy*Bv(:,i),0,N4,N4)*Av_vy*Bv;
            yhat_vy(i,:)   = B4*spdiags(Iv_vy*Bv(:,i),0,N4,N4)*yAv_vy;            
            
        end
        
        
        % evaluate convection
        convu = B1*( (yIu_ux.*(Au_ux*Bu*ru)) + (yIu_ux.*yAu_ux) ) + ...
                B2*( (yIv_uy.*(Au_uy*Bu*ru)) + (yIv_uy.*yAu_uy) );
        convv = B3*( (yIu_vx.*(Av_vx*Bv*rv)) + (yIu_vx.*yAv_vx) ) + ...
                B4*( (yIv_vy.*(Av_vy*Bv*rv)) + (yIv_vy.*yAv_vy) );
            
        for i=1:M
            convu  = convu + ru(i)*(squeeze(Chat_ux(i,:,:))*ru + yhat_ux(i,:)') + ...
                             rv(i)*(squeeze(Chat_uy(i,:,:))*ru + yhat_uy(i,:)');                     
            convv  = convv + ru(i)*(squeeze(Chat_vx(i,:,:))*rv + yhat_vx(i,:)') + ...
                             rv(i)*(squeeze(Chat_vy(i,:,:))*rv + yhat_vy(i,:)');                     
                         
        end
        
        
%         cu     = uh;
%         cv     = vh;
%         
%         u_ux   = Au_ux*uh+yAu_ux;                 % u at ux
%         uf_ux  = Iu_ux*cu+yIu_ux;                 % ubar at ux
%         du2dx  = Cux*(uf_ux.*u_ux);
%         
%         u_uy   = Au_uy*uh+yAu_uy;                 % u at uy
%         vf_uy  = Iv_uy*cv+yIv_uy;                 % vbar at uy
%         duvdy  = Cuy*(vf_uy.*u_uy);
%         
%         v_vx   = Av_vx*vh+yAv_vx;                 % v at vx
%         uf_vx  = Iu_vx*cu+yIu_vx;                 % ubar at vx
%         duvdx  = Cvx*(uf_vx.*v_vx);
%         
%         v_vy   = Av_vy*vh+yAv_vy;                 % v at vy
%         vf_vy  = Iv_vy*cv+yIv_vy;                 % vbar at vy
%         dv2dy  = Cvy*(vf_vy.*v_vy);
%         
%         convu  = du2dx + duvdy;
%         convv  = duvdx + dv2dy;
%         
%         if (getJacobian==1)
%             Newton     = options.solversettings.Newton_factor;
%             N1 = options.grid.N1;
%             N2 = options.grid.N2;
%             N3 = options.grid.N3;
%             N4 = options.grid.N4;
%             
%             C1         = Cux*spdiags(uf_ux,0,N1,N1);
%             C2         = Cux*spdiags(u_ux,0,N1,N1)*Newton;
%             Conv_ux_11 = C1*Au_ux + C2*Iu_ux;
%             
%             C1         = Cuy*spdiags(vf_uy,0,N2,N2);
%             C2         = Cuy*spdiags(u_uy,0,N2,N2)*Newton;
%             Conv_uy_11 = C1*Au_uy;
%             Conv_uy_12 = C2*Iv_uy;
%             
%             Jacu       = [Conv_ux_11 + Conv_uy_11 Conv_uy_12];
%             
%             C1         = Cvx*spdiags(uf_vx,0,N3,N3);
%             C2         = Cvx*spdiags(v_vx,0,N3,N3)*Newton;
%             Conv_vx_21 = C2*Iu_vx;
%             Conv_vx_22 = C1*Av_vx;
%             
%             C1         = Cvy*spdiags(vf_vy,0,N4,N4);
%             C2         = Cvy*spdiags(v_vy,0,N4,N4)*Newton;
%             Conv_vy_22 = C1*Av_vy + C2*Iv_vy;
%             
%             Jacv       = [Conv_vx_21 Conv_vx_22 + Conv_vy_22];
%         end
%         
        
    else
        error('not implemented');
    end
 

else
    
    error('not implemented');
    
end