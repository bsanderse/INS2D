function dvdt = dvdtBC(x,y,t,Re)

 dvdt = zeros(length(x)*length(y),1);

%     if (length(y)==1 && (abs(y-2.25)<eps || abs(y-0.25)<eps))
%         dvdt = (-2*(pi^2)/Re)*cos(pi*x)*sin(pi*y)*exp(-2*(pi^2)*t/Re);
%     else
%         dvdt = zeros(length(x)*length(y),1);
%     end
%     van Kan
%     if (t>0 & abs(y-1)<1e-12 & length(y)==1)
%         dvdt = -sin(pi*(x.^3-3*x.^2+3*x))*exp(1-1/t)/(t^2);
%     else
%         dvdt = zeros(length(x)*length(y),1);
%     end

%     if ( (length(x)==1 && abs(x-0)<1e-10)) % || (length(y)==1 && abs(y-0)<1e-10) )
%         a     = (pi/6);
%         f     = 0.5;
%         dvdt  = cos(a*sin(f*t))*a*cos(f*t)*f*ones(length(x)*length(y),1);
%     else
%         dvdt  = zeros(length(x)*length(y),1);
%     end
end