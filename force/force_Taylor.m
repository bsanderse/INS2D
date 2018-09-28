% force to obtain Taylor vortex when no pressure is used
% (force acts like pressure)

Fx = pi*exp(-4*(pi^2)*t/Re) * sin(pi*xu).*cos(pi*xu);

Fy = pi*exp(-4*(pi^2)*t/Re) * sin(pi*yv).*cos(pi*yv);


Fx = Omu.*Fx(:);
Fy = Omv.*Fy(:);