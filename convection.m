function [convu, convv] = convection(uh,vh,t,options)
% evaluate convective terms

order4     = options.discretization.order4;
regularize = options.case.regularize;

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

%%

% assumes uh, vh, cu, cv
% normally cu = uh, cv = vh

if (order4==0)
    
    %% no regularization
    if (regularize == 0)

        cu     = uh;
        cv     = vh;

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
    end


    %% leray
    if (regularize == 1)

        % filter the convecting field
        cu     = filter_convection(uh,Diffu_f,yDiffu_f,alfa); %uh + (alfa^2)*Re*(Diffu*uh + yDiffu);
        cv     = filter_convection(vh,Diffv_f,yDiffv_f,alfa);    

        % divergence of filtered velocity field; should be zero!
        maxdiv_f(n) = max(abs(M*[cu;cv]+yM));


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

    end


    %% C2
    if (regularize==2)


        % filter both convecting and convected velocity
        uh_f = filter_convection(uh,Diffu_f,yDiffu_f,alfa); %uh + (alfa^2)*Re*(Diffu*uh + yDiffu);
        vh_f = filter_convection(vh,Diffv_f,yDiffv_f,alfa);
        cu   = uh_f;
        cv   = vh_f;

        % divergence of filtered velocity field; should be zero!
        maxdiv_f(n) = max(abs(M*[uh_f;vh_f]+yM));

        u_ux   = Au_ux*uh_f+yAu_ux;                 % u at ux
        uf_ux  = Iu_ux*cu+yIu_ux;                 % ubar at ux
        du2dx  = Cux*(uf_ux.*u_ux);    

        u_uy   = Au_uy*uh_f+yAu_uy;                 % u at uy
        vf_uy  = Iv_uy*cv+yIv_uy;                 % vbar at uy
        duvdy  = Cuy*(vf_uy.*u_uy);

        v_vx   = Av_vx*vh_f+yAv_vx;                 % v at vx
        uf_vx  = Iu_vx*cu+yIu_vx;                 % ubar at vx
        duvdx  = Cvx*(uf_vx.*v_vx);

        v_vy   = Av_vy*vh_f+yAv_vy;                 % v at vy
        vf_vy  = Iv_vy*cv+yIv_vy;                 % vbar at vy    
        dv2dy  = Cvy*(vf_vy.*v_vy);

        convu  = du2dx + duvdy;
        convv  = duvdx + dv2dy;

        convu = filter_convection(convu,Diffu_f,yDiffu_f,alfa);
        convv = filter_convection(convv,Diffv_f,yDiffv_f,alfa);

    end

    %% C4
    if (regularize==4)

        % C4 consists of 3 terms:
        % C4 = conv(filter(u),filter(u)) + filter(conv(filter(u),u') +
        %      filter(conv(u',filter(u)))
        % where u' = u - filter(u)

       % filter both convecting and convected velocity
        uh_f = filter_convection(uh,Diffu_f,yDiffu_f,alfa); %uh + (alfa^2)*Re*(Diffu*uh + yDiffu);
        vh_f = filter_convection(vh,Diffv_f,yDiffv_f,alfa);
        duh  = uh - uh_f;
        dvh  = vh - vh_f;

        % divergence of filtered velocity field; should be zero!
        maxdiv_f(n) = max(abs(M*[uh_f;vh_f]+yM));

        % first term: C(filter(u),filter(u))
        u_ux   = Au_ux*uh_f+yAu_ux;                 % u at ux
        uf_ux  = Iu_ux*uh_f+yIu_ux;                 % ubar at ux
        du2dx  = Cux*(uf_ux.*u_ux);    

        u_uy   = Au_uy*uh_f+yAu_uy;                 % u at uy
        vf_uy  = Iv_uy*vh_f+yIv_uy;                 % vbar at uy
        duvdy  = Cuy*(vf_uy.*u_uy);

        v_vx   = Av_vx*vh_f+yAv_vx;                 % v at vx
        uf_vx  = Iu_vx*uh_f+yIu_vx;                 % ubar at vx
        duvdx  = Cvx*(uf_vx.*v_vx);

        v_vy   = Av_vy*vh_f+yAv_vy;                 % v at vy
        vf_vy  = Iv_vy*vh_f+yIv_vy;                 % vbar at vy    
        dv2dy  = Cvy*(vf_vy.*v_vy);

        convu1  = du2dx + duvdy;
        convv1  = duvdx + dv2dy;

        % second term: C(filter(u),u')
        u_ux   = Au_ux*duh+yAu_ux;                 % u at ux
        uf_ux  = Iu_ux*uh_f+yIu_ux;                 % ubar at ux
        du2dx  = Cux*(uf_ux.*u_ux);    

        u_uy   = Au_uy*duh+yAu_uy;                 % u at uy
        vf_uy  = Iv_uy*vh_f+yIv_uy;                 % vbar at uy
        duvdy  = Cuy*(vf_uy.*u_uy);

        v_vx   = Av_vx*dvh+yAv_vx;                 % v at vx
        uf_vx  = Iu_vx*uh_f+yIu_vx;                 % ubar at vx
        duvdx  = Cvx*(uf_vx.*v_vx);

        v_vy   = Av_vy*dvh+yAv_vy;                 % v at vy
        vf_vy  = Iv_vy*vh_f+yIv_vy;                 % vbar at vy    
        dv2dy  = Cvy*(vf_vy.*v_vy);

        convu2  = du2dx + duvdy;
        convv2  = duvdx + dv2dy;

        % third term: C(u',filter(u))
        u_ux   = Au_ux*uh_f+yAu_ux;                 % u at ux
        uf_ux  = Iu_ux*duh+yIu_ux;                 % ubar at ux
        du2dx  = Cux*(uf_ux.*u_ux);    

        u_uy   = Au_uy*uh_f+yAu_uy;                 % u at uy
        vf_uy  = Iv_uy*dvh+yIv_uy;                 % vbar at uy
        duvdy  = Cuy*(vf_uy.*u_uy);

        v_vx   = Av_vx*vh_f+yAv_vx;                 % v at vx
        uf_vx  = Iu_vx*duh+yIu_vx;                 % ubar at vx
        duvdx  = Cvx*(uf_vx.*v_vx);

        v_vy   = Av_vy*vh_f+yAv_vy;                 % v at vy
        vf_vy  = Iv_vy*dvh+yIv_vy;                 % vbar at vy    
        dv2dy  = Cvy*(vf_vy.*v_vy);

        convu3  = du2dx + duvdy;
        convv3  = duvdx + dv2dy;    

        convu = convu1 + filter_convection(convu2+convu3,Diffu_f,yDiffu_f,alfa);
        convv = convv1 + filter_convection(convv2+convv3,Diffv_f,yDiffv_f,alfa);


    end

elseif (order4==1)

    alfa   = options.discretization.alfa;
    
    Au_ux3 = options.discretization.Au_ux3;
    Au_uy3 = options.discretization.Au_uy3;
    Av_vx3 = options.discretization.Av_vx3;
    Av_vy3 = options.discretization.Av_vy3;
    yAu_ux3 = options.discretization.yAu_ux3;
    yAu_uy3 = options.discretization.yAu_uy3;
    yAv_vx3 = options.discretization.yAv_vx3;
    yAv_vy3 = options.discretization.yAv_vy3;    
    
    Iu_ux3 = options.discretization.Iu_ux3;
    Iv_uy3 = options.discretization.Iv_uy3;
    Iu_vx3 = options.discretization.Iu_vx3;
    Iv_vy3 = options.discretization.Iv_vy3;
    yIu_ux3 = options.discretization.yIu_ux3;
    yIv_uy3 = options.discretization.yIv_uy3;
    yIu_vx3 = options.discretization.yIu_vx3;
    yIv_vy3 = options.discretization.yIv_vy3;   
    
    Cux3 = options.discretization.Cux3;
    Cuy3 = options.discretization.Cuy3;
    Cvx3 = options.discretization.Cvx3;
    Cvy3 = options.discretization.Cvy3;
    
    cu     = uh;
    cv     = vh;

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

    u_ux3  = Au_ux3*uh+yAu_ux3;                 % u at ux
    uf_ux3 = Iu_ux3*cu+yIu_ux3;                 % ubar at ux
    du2dx3 = Cux3*(uf_ux3.*u_ux3);    

    u_uy3  = Au_uy3*uh+yAu_uy3;                 % u at uy
    vf_uy3 = Iv_uy3*cv+yIv_uy3;                 % vbar at uy
    duvdy3 = Cuy3*(vf_uy3.*u_uy3);

    v_vx3  = Av_vx3*vh+yAv_vx3;                 % v at vx
    uf_vx3 = Iu_vx3*cu+yIu_vx3;                 % ubar at vx
    duvdx3 = Cvx3*(uf_vx3.*v_vx3);

    v_vy3  = Av_vy3*vh+yAv_vy3;                 % v at vy
    vf_vy3 = Iv_vy3*cv+yIv_vy3;                 % vbar at vy    
    dv2dy3 = Cvy3*(vf_vy3.*v_vy3);

    
    convu = alfa*du2dx - du2dx3 + alfa*duvdy - duvdy3;
    convv = alfa*duvdx - duvdx3 + alfa*dv2dy - dv2dy3;
    
end
    
    
end