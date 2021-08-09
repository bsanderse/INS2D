clear all; close all;

% f = @(y,t) sin(y-t); %requires two basis modes for exact representation
% f = @(y,t) sin(y.^2-t); %2 modes
% f = @(y,t) sin((y-t).^2); %requires approx 50 modes
f = @(y,t) cos(pi/6*sin(y-t/2)); %thisisit

dx = 100;
dt = 200;
y = linspace(-2,2,dx);
ts = linspace(0,4*pi,dt);

X = zeros(dx,dt);
figure
for j = 1:numel(ts)
    t = ts(j);
    f_val = f(y,t);
    X(:,j) = f_val';
    plot3(y,f(y,t),t*ones(dx,1))
%     plot(y,f(y,t))
    hold on
end

xlabel('y')
ylabel('f')
zlabel('z')

[U,S,V] = svd(X);
figure; surf(U)
figure; loglog(diag(S))
Ms = [2 4 10 20 40];
for k = 1:length(Ms)
    M = Ms(k)
    phi = U(:,1:M);
    a = phi'*X;
    figure; surf(phi*a)
    norm(X-phi*a)
end