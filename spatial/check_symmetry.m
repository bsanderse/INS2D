function [flag, symmetry_error] = check_symmetry(uh,vh,t,options)
% check symmetry of convection operator
% flag = 0: no symmetry error
% flag = 1: symmetry error

eps = 1e-14;
fcw = options.output.fcw;

BC  = options.BC;

Cux = options.discretization.Cux;
Cuy = options.discretization.Cuy;
Cvx = options.discretization.Cvx;
Cvy = options.discretization.Cvy;

Au_ux = options.discretization.Au_ux;
Au_uy = options.discretization.Au_uy;
Av_vx = options.discretization.Av_vx;
Av_vy = options.discretization.Av_vy;

Iu_ux = options.discretization.Iu_ux;
Iv_uy = options.discretization.Iv_uy;
Iu_vx = options.discretization.Iu_vx;
Iv_vy = options.discretization.Iv_vy;

yIu_ux = options.discretization.yIu_ux;
yIv_uy = options.discretization.yIv_uy;
yIu_vx = options.discretization.yIu_vx;
yIv_vy = options.discretization.yIv_vy;

N1 = options.grid.N1;
N2 = options.grid.N2;
N3 = options.grid.N3;
N4 = options.grid.N4;

Cu = Cux*spdiags(Iu_ux*uh+yIu_ux,0,N1,N1)*Au_ux + ...
     Cuy*spdiags(Iv_uy*vh+yIv_uy,0,N2,N2)*Au_uy;
Cv = Cvx*spdiags(Iu_vx*uh+yIu_vx,0,N3,N3)*Av_vx + ...
     Cvy*spdiags(Iv_vy*vh+yIv_vy,0,N4,N4)*Av_vy;

error_u = max2d(abs(Cu+Cu'));
error_v = max2d(abs(Cv+Cv'));

symmetry_error = max(error_u,error_v);

flag = 0;
if (symmetry_error > eps)
    if (~strcmp(BC.u.right,'pres') && ~strcmp(BC.u.left,'pres'))
        fprintf(fcw,[num2str(error_u) '\n']);
        flag = 1;
    end
    if (~strcmp(BC.v.low,'pres') && ~strcmp(BC.v.up,'pres'))
        fprintf(fcw,[num2str(error_v) '\n']);
        flag = 1;
    end
end

        
%     if (max2d(abs(Cu+Cu'))>1e-12 && ...
%           ~strcmp(BC.u.right,'pres') && ~strcmp(BC.u.left,'pres') )
%         disp('warning: convection operator u not skew-symmetric');
%     end
%     if (max2d(abs(Cv+Cv'))>1e-12  && ...
%           ~strcmp(BC.v.low,'pres') && ~strcmp(BC.v.up,'pres') )
%         disp('warning: convection operator v not skew-symmetric');
%     end

end