function u = uBC(x,y,t,Re)

global x1 x2 y1 y2;


u = zeros(length(x)*length(y),1);

%inflow
% if ( (length(x)==1 && abs(x-x1)<1e-10) || (length(y)==1 && abs(y-y1)<1e-10) )
% %     f     = 0.5;
% %     theta = (pi/6)*sin(f*t);
% %     u     = cos(theta)*ones(length(x)*length(y),1);
% %     u = ones(length(y),1);
%     u = cosd(30)*ones(length(x)*length(y),1);
% %     u = (2 + cos(t))*ones(length(y),1);
% 
% %     u = u + (y<=1);
% % elseif ( (length(y)==1 && abs(y-0)<1e-10) )
% %     u = ones(length(x),1);
% %     u = u + (x<=1);
% else
%     u = zeros(length(x)*length(y),1);
% end

% AD
% if (length(x)==1 && abs(x-x1)<1e-10)
% %     Ct = 0.001;
% %     u = 1 + 0.5*Ct*1 / (2*pi) * (atan((1/2-y)/x) + atan((1/2+y)/x)) ;
%     u = ones(length(y),1);
% else
%     u = zeros(length(x)*length(y),1);
% end
%  

% manufactured stuff
%     if (length(y)==1)
%         y=y*ones(length(x),1);
%     end
%     if (length(x)==1)
%         x=x*ones(length(y),1);
%     end
%     u = -(137/36)*x+2*y+x.*y+(17/8)*x.^2+(17/8)*y.^2 ...
%               -(26/9)*x.^3-2*y.^3+(13/9)*x.^4+y.^4;

% Taylor:
%       u = -sin(pi*x)*cos(pi*y)*exp(-2*(pi^2)*t/Re);
%     u = -(pi^2)*x*cos(2*pi*y);
%     if (length(y)==1 && (abs(y-2.25)<eps || abs(y-0.25)<eps))
%         u = -sin(pi*x)*cos(pi*y)*exp(-2*(pi^2)*t/Re);
%         u = -sin(pi*x)*ones(length(y),1); %*y.^2;
%      if (length(x)==1)        
%          u = -pi*cos(pi*x).*cos(2*pi*y);
%          u = zeros(length(y),1);
%      else
%          u = 2*pi*sin(pi*x).*sin(2*pi*y);
%          u = pi*sin(pi*y)*ones(length(x),1);
%         u = -cos(pi*y)*ones(length(x),1);       
%         u = -sin(pi*x)*cos(2*pi*y);
%      end
%     else
%         u = zeros(length(x)*length(y),1);
%     end

% LDC Shih
%     if (length(y)==1 && abs(y-1)<1e-10)
% %         u = 16*(x.^4-2*x.^3+x.^2);
%         u = -ones(length(x)*length(y),1);
% %     elseif (length(y)==1 && abs(y-0)<1e-10)
% %         u = -16*(x.^4-2*x.^3+x.^2);        
% %         u = -ones(length(x)*length(y),1);
%     else
%         u = zeros(length(x)*length(y),1);
%     end

% Kovasznay
%     l = Re/2-sqrt(Re^2/4+4*pi^2);
% if ( (length(x)==1 && abs(x-0.1)<1e-10) || ...
%      (length(y)==1 && abs(y-0.4)<1e-10) || ...
%      (length(y)==1 && abs(y+0.4)<1e-10) )
%     u = 1-exp(l*x)*cos(2*pi*y);
% else
%     u = zeros(length(x)*length(y),1);
% end


% BFS
%       if (length(x)==1 && x==0)
%          u = (y>=0).*(24*y.*(1/2-y)); 
%       end
    
% if (abs(x)<=1e-10)
%     u =(0.5*(-tanh((y-0.75)/0.02)+tanh((y-0.25)/0.02))+1)*(1-exp(-(t)^2));
%     u = ones(length(x)*length(y),1);%*(1-exp(-(t)^2));
%     u = sin(2*pi*y)*sin(t)+2;
% else
%     u = zeros(length(x)*length(y),1);    
% end

% LDC
    if (length(y)==1 && abs(y-1)<1e-10)
        u = ones(length(x)*length(y),1);
    end

% van Kan
%     u = zeros(length(x)*length(y),1);    
    
% unsteady wake
% if (x==0 & length(x)==1)
%     u = (1.5*(y>0)+0.5)*( ...
%           exp(-100*(t-0.5)^2)*(t<0.5) + 1*(t>=0.5));
%           exp(-100*(t-1.5)^2)*(t<1.5) + );
%     u = (y>0)*(cos(t)+2) + (y<=0)*(2-cos(t));
% else
%     u = zeros(length(x)*length(y),1); 
% end