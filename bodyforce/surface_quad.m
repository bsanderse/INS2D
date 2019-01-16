eps = 0.000;
d   = 1e-14;

U   = [eps;1;0];
dt  = 1;

x11=[0;0;1];
x21=[0;0.5;1];

x12=x11+U*dt;
x22=x21+U*dt;

x=[x11(1) x21(1) x22(1) x12(1)];
y=[x11(2) x21(2) x22(2) x12(2)];

fill(x,y,'r');
axis([0 1 0 2]);

% x11=[0;0;1];
% x21=[0;1;1];
% 
% x22=[1;2;1];
% x12=[1;1;1];

ez = [0;0;1];

xa = 0.5*(x11+x21);
xb = 0.5*(x21+x22);
xc = 0.5*(x22+x12);
xd = 0.5*(x12+x11);

na = cross(x11-x21,ez);
na = na/norm(na);
nb = cross(x21-x22,ez);
nb = nb/norm(nb);
nc = cross(x22-x12,ez);
nc = nc/norm(nc);
nd = cross(x12-x11,ez);
nd = nd/norm(nd);

% V = dot(xa,na)*norm(x21-x11) + dot(xb,nb)*norm(x22-x21) + ...
%     dot(xc,nc)*norm(x12-x22) + dot(xd,nd)*norm(x11-x12)

% as cross product:
V = norm(cross(x21-x11,x12-x11+d*na))
% as determinant:
% V = det( [x11';x12';x21'] ) + d

un = abs(dot(x12-x11+d*na,na))/dt

V/un