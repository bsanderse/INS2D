function v = vBC(x,y,t,Re)

global x1 x2 y1 y2;

v = zeros(length(x)*length(y),1);

% rotated LDC
%     if (length(x)==1 && abs(x-1)<1e-10)
% %         u = 16*(x.^4-2*x.^3+x.^2);
%         v = ones(length(x)*length(y),1);
%     else
%         v = zeros(length(x)*length(y),1);
%     end
% 
% if (length(x)==1 && abs(x-1)<1e-10)
%     v = -ones(length(x)*length(y),1);
% elseif (length(x)==1 && abs(x-0)<1e-10)
%     v = ones(length(x)*length(y),1);
% else
%     v = zeros(length(x)*length(y),1);
% end


% if ( (length(x)==1 && abs(x-0)<1e-10)) % || (length(y)==1 && abs(y-0)<1e-10) )
% %     f     = 0.5;
% %     theta = (pi/6)*sin(f*t);
% %     u     = cos(theta)*ones(length(x)*length(y),1);
%     v = ones(length(y),1);
% %     u(ceil(length(u)/2)) = 2;
% %     u = u + (y<=5.5 & y>=4.5);
%     v = v + (y<=1);
% %     keyboard
% end
% if ( (length(y)==1 && abs(y-0)<1e-10) )
%     v = ones(length(x),1);
%     v = v + (x<=1);
% end

% Taylor
%  v = cos(pi*x)*sin(pi*y)*exp(-2*(pi^2)*t/Re);
%    v = zeros(length(x)*length(y),1);
%     if (length(y)==1 && (abs(y-2.25)<eps || abs(y-0.25)<eps))
%         v = cos(pi*x)*sin(pi*y)*exp(-2*(pi^2)*t/Re);
%     else
%         v = zeros(length(x)*length(y),1);        
%     end
        
% if ( (length(x)==1 && abs(x-x1)<1e-10) || (length(y)==1 && abs(y-y1)<1e-10) )
%     f     = 0.5;
%     theta = (pi/6)*sin(f*t);
%     v     = sin(theta)*ones(length(x)*length(y),1);
%     v = sind(30)*ones(length(x)*length(y),1);
% else
%     v = zeros(length(x)*length(y),1);
% end
% % van Kan:
%     if (t>0 & abs(y-1)<1e-12 & length(y)==1)
%         v = -sin(pi*(x.^3-3*x.^2+3*x))*exp(1-1/t);
%     else
%         v = zeros(length(x)*length(y),1);
%     end

% steady van Kan
%     if (abs(y-1)<1e-12 & length(y)==1)
% %         v = -sin(pi*(x.^3-3*x.^2+3*x));
% %         v = -(27/4)*x.*(x-1).^2;
% %           v = -(256/27)*(-x.^4+3*x.^3-3*x.^2+x);
%           xm = 1-2^(-1/3);
%           v  = -x.*(x-1).^4/(xm*(xm-1)^4);
%     else
%         v = zeros(length(x)*length(y),1);
%     end        

% Poiseuille like
%     if (length(y)==1)
% %         v = 0.2*(sin(pi*x).^2-1/2);
%         v = 0.2*sin(2*pi*x)*cos(pi*y);
%     else
%         v = zeros(length(x)*length(y),1);
%     end        
    
% if (abs(y)<=1e-10)
%     v = ones(length(x)*length(y),1);
% else
%     v = zeros(length(x)*length(y),1);    
% end
% Kovasznay
% l = Re/2-sqrt(Re^2/4+4*pi^2);
% if ( (length(x)==1 && abs(x-0.1)<1e-10) || ...
%         (length(y)==1 && abs(y-0.4)<1e-10) || ...
%         (length(y)==1 && abs(y+0.4)<1e-10) )
%     v = l/(2*pi) * exp(l*x)*sin(2*pi*y);
% else
%     % dvdx is nonzero at outlet
%     v = (l^2)/(2*pi)*exp(l*x)*sin(2*pi*y);
% %     v = zeros(length(x)*length(y),1);
% end

end