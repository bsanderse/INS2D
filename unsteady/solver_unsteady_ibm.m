x_ibm = Nvx_in;
y_ibm = Nvy_in/2;
n_ibm = floor(1/hx(1))-1;
n_ibm_start = y_ibm*Nvx_in + Nvx_in/4;
% ibm according to Colonius

Ex = spalloc(n_ibm,Nu,0);

E1 = spalloc(n_ibm,n_ibm_start,0);
E2 = spdiags(ones(n_ibm,1),0,n_ibm,n_ibm);
E3 = spalloc(n_ibm,Nv-n_ibm_start-n_ibm,0);
Ey = [E1 E2 E3];

Mx = [Mx; Ex];
My = [My; Ey];
M  = [Mx My];
Gx = [Gx Ex'];
Gy = [Gy Ey'];
G  = [Gx; Gy];

Np = Np + n_ibm;
