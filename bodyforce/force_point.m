% point force
xmid       = find(x>=L_x/2,1);
ymid       = floor(Nvy_in/2)+1;
Gamma      = -0.11;
y1Dx       = zeros(Nvx_in,1);
y1Dx(xmid) = 1;
y1Dy       = zeros(Nvy_in,1);
y1Dy(ymid) = Gamma;

Fy         = kron(y1Dy,y1Dx);
