function dudt = dudtBC(x,y,t,Re)

 dudt = zeros(length(x)*length(y),1);

% Taylor:
%     if (length(y)==1 && (abs(y-2.25)<eps || abs(y-0.25)<eps))
%         dudt = (2*(pi^2)/Re)*sin(pi*x)*cos(pi*y)*exp(-2*(pi^2)*t/Re);
%     else
%         dudt = zeros(length(x)*length(y),1);
%     end
%    if (x==0)
%     dudt = sin(2*pi*y)*cos(t);
%     dudt =(0.5*(-tanh((y-0.75)/0.02)+tanh((y-0.25)/0.02))+1)*exp(-(t)^2)*2*(t);
%     dudt = ones(length(x)*length(y),1)*exp(-(t)^2)*2*(t);
%    else
%     dudt = zeros(length(x)*length(y),1);
%    end

% unsteady wake
% if (abs(x-0)<1e-10 & length(x)==1)
% 
%     a     = (pi/6);
%     f     = 0.5;
%     dudt  = -sin(a*sin(f*t))*a*cos(f*t)*f*ones(length(x)*length(y),1);
% 
% else
%     dudt = zeros(length(x)*length(y),1); 
% end

end