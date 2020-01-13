function [convu, convv, Jacu, Jacv] = convection(V,C,t,options,getJacobian)
% evaluate convective terms and, optionally, Jacobians
% V: velocity field
% C: 'convection' field: e.g. d(c_x u)/dx + d(c_y u)/dy; usually c_x = u,
% c_y=v

alfa     = options.discretization.alfa;
order4   = options.discretization.order4;

regularize = options.case.regularize;

Nu = options.grid.Nu;
Nv = options.grid.Nv;
NV = options.grid.NV;

indu = options.grid.indu;
indv = options.grid.indv;

Jacu = spalloc(Nu,NV,0);
Jacv = spalloc(Nv,NV,0);

%%
uh = V(indu);
vh = V(indv);

cu = C(indu);
cv = C(indv);


switch regularize

    case 0


        % no regularization
        [convu,convv,Jacu,Jacv] = convection_components(C,V,options,getJacobian,0);

        if (order4==1)
            [convu3,convv3,Jacu3,Jacv3] = convection_components(C,V,options,getJacobian,1);
            convu = alfa*convu - convu3;
            convv = alfa*convv - convv3;
            Jacu  = alfa*Jacu - Jacu3;
            Jacv  = alfa*Jacv - Jacv3;
        end

    case 1
        % Leray
        % TODO: needs finishing

        % filter the convecting field
        cu_f     = filter_convection(cu,Diffu_f,yDiffu_f,alfa); %uh + (alfa^2)*Re*(Diffu*uh + yDiffu);
        cv_f     = filter_convection(cv,Diffv_f,yDiffv_f,alfa);

        C_filtered = [cu_f;cv_f];

        % divergence of filtered velocity field; should be zero!
        maxdiv_f = max(abs(M*C_filtered+yM));

        [convu,convv,Jacu,Jacv] = convection_components(C_filtered,V,options,getJacobian);

    case 2
        %% C2

        cu_f     = filter_convection(cu,Diffu_f,yDiffu_f,alfa); %uh + (alfa^2)*Re*(Diffu*uh + yDiffu);
        cv_f     = filter_convection(cv,Diffv_f,yDiffv_f,alfa);

        uh_f     = filter_convection(uh,Diffu_f,yDiffu_f,alfa); %uh + (alfa^2)*Re*(Diffu*uh + yDiffu);
        vh_f     = filter_convection(vh,Diffv_f,yDiffv_f,alfa);

        C_filtered = [cu_f;cv_f];
        V_filtered = [uh_f;vh_f];

        % divergence of filtered velocity field; should be zero!
        maxdiv_f = max(abs(M*C_filtered+yM));

        [convu,convv,Jacu,Jacv] = convection_components(C_filtered,V_filtered,options,getJacobian);

        convu = filter_convection(convu,Diffu_f,yDiffu_f,alfa);
        convv = filter_convection(convv,Diffv_f,yDiffv_f,alfa);


    case 4

        % C4 consists of 3 terms:
        % C4 = conv(filter(u),filter(u)) + filter(conv(filter(u),u') +
        %      filter(conv(u',filter(u)))
        % where u' = u - filter(u)

        % filter both convecting and convected velocity
        uh_f = filter_convection(uh,Diffu_f,yDiffu_f,alfa); %uh + (alfa^2)*Re*(Diffu*uh + yDiffu);
        vh_f = filter_convection(vh,Diffv_f,yDiffv_f,alfa);

        V_filtered = [uh_f;vh_f];

        dV   = V - V_filtered;

        cu_f = filter_convection(cu,Diffu_f,yDiffu_f,alfa); %uh + (alfa^2)*Re*(Diffu*uh + yDiffu);
        cv_f = filter_convection(cv,Diffv_f,yDiffv_f,alfa);

        C_filtered = [cu_f;cv_f];
        dC   = C - C_filtered;


        % divergence of filtered velocity field; should be zero!
        maxdiv_f(n) = max(abs(M*V_filtered+yM));

        [convu1,convv1,Jacu,Jacv] = convection_components(C_filtered,V_filtered,options,getJacobian);

        [convu2,convv2,Jacu,Jacv] = convection_components(C_filtered,dV,options,getJacobian);

        [convu3,convv3,Jacu,Jacv] = convection_components(dC,V_filtered,options,getJacobian);

        convu = convu1 + filter_convection(convu2+convu3,Diffu_f,yDiffu_f,alfa);
        convv = convv1 + filter_convection(convv2+convv3,Diffv_f,yDiffv_f,alfa);


end
    



end

function [convu,convv,Jacu,Jacv] = convection_components(C,V,options,getJacobian,order4)

if (nargin<=4)
    order4 = 0;
end

if (order4 == 0)
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
elseif (order4 == 1)
    Cux = options.discretization.Cux3;
    Cuy = options.discretization.Cuy3;
    Cvx = options.discretization.Cvx3;
    Cvy = options.discretization.Cvy3;

    Au_ux = options.discretization.Au_ux3;
    Au_uy = options.discretization.Au_uy3;
    Av_vx = options.discretization.Av_vx3;
    Av_vy = options.discretization.Av_vy3;

    yAu_ux = options.discretization.yAu_ux3;
    yAu_uy = options.discretization.yAu_uy3;
    yAv_vx = options.discretization.yAv_vx3;
    yAv_vy = options.discretization.yAv_vy3;

    Iu_ux = options.discretization.Iu_ux3;
    Iv_uy = options.discretization.Iv_uy3;
    Iu_vx = options.discretization.Iu_vx3;
    Iv_vy = options.discretization.Iv_vy3;

    yIu_ux = options.discretization.yIu_ux3;
    yIv_uy = options.discretization.yIv_uy3;
    yIu_vx = options.discretization.yIu_vx3;
    yIv_vy = options.discretization.yIv_vy3;   
end

indu = options.grid.indu;
indv = options.grid.indv;

Nu = options.grid.Nu;
Nv = options.grid.Nv;
NV = options.grid.NV;

Jacu = spalloc(Nu,NV,0);
Jacv = spalloc(Nv,NV,0);

uh = V(indu);
vh = V(indv);

cu = C(indu);
cv = C(indv);

u_ux   = Au_ux*uh+yAu_ux;                 % u at ux
uf_ux  = Iu_ux*cu+yIu_ux;                 % ubar at ux
du2dx  = Cux*(uf_ux.*u_ux);

u_uy   = Au_uy*uh+yAu_uy;                 % u at uy
vf_uy  = Iv_uy*cv+yIv_uy;                 % vbar at uy
duvdy  = Cuy*(vf_uy.*u_uy);

v_vx   = Av_vx*vh+yAv_vx;                 % v at vx
uf_vx  = Iu_vx*cu+yIu_vx;                 % ubar at vx
duvdx  = Cvx*(uf_vx.*v_vx);

v_vy   = Av_vy*vh+yAv_vy;                 % v at vy
vf_vy  = Iv_vy*cv+yIv_vy;                 % vbar at vy
dv2dy  = Cvy*(vf_vy.*v_vy);

convu  = du2dx + duvdy;
convv  = duvdx + dv2dy;

if (getJacobian==1)
    Newton     = options.solversettings.Newton_factor;
    N1 = length(u_ux); %options.grid.N1;
    N2 = length(u_uy); %options.grid.N2;
    N3 = length(v_vx); %options.grid.N3;
    N4 = length(v_vy); %options.grid.N4;
    
    
    %% convective terms, u-component
    % c^n * u^(n+1), c=u
    C1         = Cux*spdiags(uf_ux,0,N1,N1);
    C2         = Cux*spdiags(u_ux,0,N1,N1)*Newton;
    Conv_ux_11 = C1*Au_ux + C2*Iu_ux;
    
    C1         = Cuy*spdiags(vf_uy,0,N2,N2);
    C2         = Cuy*spdiags(u_uy,0,N2,N2)*Newton;
    Conv_uy_11 = C1*Au_uy;
    Conv_uy_12 = C2*Iv_uy;
    
    Jacu       = [Conv_ux_11 + Conv_uy_11 Conv_uy_12];
    
    %% convective terms, v-component
    C1         = Cvx*spdiags(uf_vx,0,N3,N3);
    C2         = Cvx*spdiags(v_vx,0,N3,N3)*Newton;
    Conv_vx_21 = C2*Iu_vx;
    Conv_vx_22 = C1*Av_vx;
    
    C1         = Cvy*spdiags(vf_vy,0,N4,N4);
    C2         = Cvy*spdiags(v_vy,0,N4,N4)*Newton;
    Conv_vy_22 = C1*Av_vy + C2*Iv_vy;
    
    Jacv       = [Conv_vx_21 Conv_vx_22 + Conv_vy_22];
    
end

end

% old 4th order implementation:
% elseif (order4==1)
%     
%     alfa   = options.discretization.alfa;
%     
%     Au_ux3 = options.discretization.Au_ux3;
%     Au_uy3 = options.discretization.Au_uy3;
%     Av_vx3 = options.discretization.Av_vx3;
%     Av_vy3 = options.discretization.Av_vy3;
%     yAu_ux3 = options.discretization.yAu_ux3;
%     yAu_uy3 = options.discretization.yAu_uy3;
%     yAv_vx3 = options.discretization.yAv_vx3;
%     yAv_vy3 = options.discretization.yAv_vy3;
%     
%     Iu_ux3 = options.discretization.Iu_ux3;
%     Iv_uy3 = options.discretization.Iv_uy3;
%     Iu_vx3 = options.discretization.Iu_vx3;
%     Iv_vy3 = options.discretization.Iv_vy3;
%     yIu_ux3 = options.discretization.yIu_ux3;
%     yIv_uy3 = options.discretization.yIv_uy3;
%     yIu_vx3 = options.discretization.yIu_vx3;
%     yIv_vy3 = options.discretization.yIv_vy3;
%     
%     Cux3 = options.discretization.Cux3;
%     Cuy3 = options.discretization.Cuy3;
%     Cvx3 = options.discretization.Cvx3;
%     Cvy3 = options.discretization.Cvy3;
%     
%     %     cu     = uh;
%     %     cv     = vh;
%     
%     u_ux   = Au_ux*uh+yAu_ux;                 % u at ux
%     uf_ux  = Iu_ux*cu+yIu_ux;                 % ubar at ux
%     du2dx  = Cux*(uf_ux.*u_ux);
%     
%     u_uy   = Au_uy*uh+yAu_uy;                 % u at uy
%     vf_uy  = Iv_uy*cv+yIv_uy;                 % vbar at uy
%     duvdy  = Cuy*(vf_uy.*u_uy);
%     
%     v_vx   = Av_vx*vh+yAv_vx;                 % v at vx
%     uf_vx  = Iu_vx*cu+yIu_vx;                 % ubar at vx
%     duvdx  = Cvx*(uf_vx.*v_vx);
%     
%     v_vy   = Av_vy*vh+yAv_vy;                 % v at vy
%     vf_vy  = Iv_vy*cv+yIv_vy;                 % vbar at vy
%     dv2dy  = Cvy*(vf_vy.*v_vy);
%     
%     u_ux3  = Au_ux3*uh+yAu_ux3;                 % u at ux
%     uf_ux3 = Iu_ux3*cu+yIu_ux3;                 % ubar at ux
%     du2dx3 = Cux3*(uf_ux3.*u_ux3);
%     
%     u_uy3  = Au_uy3*uh+yAu_uy3;                 % u at uy
%     vf_uy3 = Iv_uy3*cv+yIv_uy3;                 % vbar at uy
%     duvdy3 = Cuy3*(vf_uy3.*u_uy3);
%     
%     v_vx3  = Av_vx3*vh+yAv_vx3;                 % v at vx
%     uf_vx3 = Iu_vx3*cu+yIu_vx3;                 % ubar at vx
%     duvdx3 = Cvx3*(uf_vx3.*v_vx3);
%     
%     v_vy3  = Av_vy3*vh+yAv_vy3;                 % v at vy
%     vf_vy3 = Iv_vy3*cv+yIv_vy3;                 % vbar at vy
%     dv2dy3 = Cvy3*(vf_vy3.*v_vy3);
%     
%     
%     convu = alfa*du2dx - du2dx3 + alfa*duvdy - duvdy3;
%     convv = alfa*duvdx - duvdx3 + alfa*dv2dy - dv2dy3;
%     
%     if (getJacobian==1)
%         Newton     = options.solversettings.Newton_factor;
%         N1 = options.grid.N1;
%         N2 = options.grid.N2;
%         N3 = options.grid.N3;
%         N4 = options.grid.N4;
%         
%         %% 2nd order operations
%         C1         = Cux*spdiags(uf_ux,0,N1,N1);
%         C2         = Cux*spdiags(u_ux,0,N1,N1)*Newton;
%         Conv_ux_11 = C1*Au_ux + C2*Iu_ux;
%         C1         = Cuy*spdiags(vf_uy,0,N2,N2);
%         C2         = Cuy*spdiags(u_uy,0,N2,N2)*Newton;
%         Conv_uy_11 = C1*Au_uy;
%         Conv_uy_12 = C2*Iv_uy;
%         Jacu1      = [Conv_ux_11 + Conv_uy_11 Conv_uy_12];
%         
%         C1         = Cvx*spdiags(uf_vx,0,N3,N3);
%         C2         = Cvx*spdiags(v_vx,0,N3,N3)*Newton;
%         Conv_vx_21 = C2*Iu_vx;
%         Conv_vx_22 = C1*Av_vx;
%         C1         = Cvy*spdiags(vf_vy,0,N4,N4);
%         C2         = Cvy*spdiags(v_vy,0,N4,N4)*Newton;
%         Conv_vy_22 = C1*Av_vy + C2*Iv_vy;
%         Jacv1       = [Conv_vx_21 Conv_vx_22 + Conv_vy_22];
%         
%         %% 4th order operations
%         C1     = Cux3*spdiags(uf_ux3,0,length(uf_ux3),length(uf_ux3));
%         C2     = Cux3*spdiags(u_ux3,0,length(u_ux3),length(u_ux3))*Newton;
%         Conv_ux_11_3 = C1*Au_ux3 + C2*Iu_ux3;
%         C1         = Cuy3*spdiags(vf_uy3,0,length(vf_uy3),length(vf_uy3));
%         C2         = Cuy3*spdiags(u_uy3,0,length(vf_uy3),length(vf_uy3))*Newton;
%         Conv_uy_11_3 = C1*Au_uy3;
%         Conv_uy_12_3 = C2*Iv_uy3;
%         Jacu3      = [Conv_ux_11_3 + Conv_uy_11_3 Conv_uy_12_3];
%         
%         
%         %%  v-component
%         C1         = Cvx3*spdiags(uf_vx3,0,length(uf_vx3),length(uf_vx3));
%         C2         = Cvx3*spdiags(v_vx3,0,length(v_vx3),length(v_vx3))*Newton;
%         Conv_vx_21_3 = C2*Iu_vx3;
%         Conv_vx_22_3 = C1*Av_vx3;
%         C1         = Cvy3*spdiags(vf_vy3,0,length(vf_vy3),length(vf_vy3));
%         C2         = Cvy3*spdiags(v_vy3,0,length(v_vy3),length(v_vy3))*Newton;
%         Conv_vy_22_3 = C1*Av_vy3 + C2*Iv_vy3;
%         Jacv3      =  [Conv_vx_21_3 Conv_vx_22_3 + Conv_vy_22_3];
%         
%         % linear combination to get the 4th order operator
%         Jacu  = alfa*Jacu1 - Jacu3;
%         Jacv  = alfa*Jacv1 - Jacv3;
%     end
%     
%     
% end