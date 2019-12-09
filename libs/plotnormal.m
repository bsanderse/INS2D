function plotnormal(x_k,y_k,vect,vecn)


plot(x_k,y_k,'kx-')

hold on

S    = sqrt((x_k(1:end-1)-x_k(2:end)).^2 + (y_k(1:end-1)-y_k(2:end)).^2);

xcol = 0.5*(x_k(1:end-1)+x_k(2:end));
ycol = 0.5*(y_k(1:end-1)+y_k(2:end));

plot(xcol,ycol,'bo');

% tangential vectors
xt = [xcol'; (xcol+0.5*S.*vect(:,1))'];
yt = [ycol'; (ycol+0.5*S.*vect(:,2))'];
plot(xt,yt,'rd-');

% normal vectors
xn = [xcol'; (xcol+0.5*S.*vecn(:,1))'];
yn = [ycol'; (ycol+0.5*S.*vecn(:,2))'];
plot(xn,yn,'gs-');
