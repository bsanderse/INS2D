% mesh comparison

Nx = 64;

xi = nonuniform_grid(1/Nx,0,1,1);


x_cos = (1/2) * (1 - cos(pi*xi));
x_exp = nonuniform_grid2(0,1/2,(2/pi)*Nx,Nx/2);
x_exp = [x_exp(1:end-1); flipud(1-x_exp)];

hx_cos = diff(x_cos);
hx_exp = diff(x_exp);

xm_cos = 0.5*(x_cos(1:end-1)+x_cos(2:end));
xm_exp = 0.5*(x_exp(1:end-1)+x_exp(2:end));

% should be approximately equal:
max(hx_cos)/min(hx_cos)
max(hx_exp)/min(hx_exp)

min(hx_cos)
0.25*(pi*(1/Nx))^2
max(hx_cos)
0.5*pi*(1/Nx)

figure
plot(xm_cos,hx_cos,'rx-')
hold on
plot(xm_exp,hx_exp,'bx-')
legend('cosine','exponential')
grid
layout(1,16,'x','\Delta x')
