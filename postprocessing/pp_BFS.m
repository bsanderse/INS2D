% post-processing
u    = reshape(uh,Nux_in,Nuy_in);
v    = reshape(vh,Nvx_in,Nvy_in);
pres = reshape(p,Npx,Npy);


% calculate length of reattachment
figure(1)

start = 8; % prevent first points to be taken
plot(x(2:end),u(:,1),'b');
hold on
plot(x(2:end),u(:,end),'r');
grid
legend(['u at y=' num2str(yu(1,1))],['u at y=' num2str(yu(1,end))])
xlabel('x'); ylabel('u')

ubottom = u(start:end,1);

% find first point that has positive velocity
il = find(min(ubottom,0)==0,1);
if (~isempty(il))
    % interpolate to obtain intersection point
    X1(j) = x(il+start-1) + (-ubottom(il-1)/(ubottom(il)-ubottom(il-1)))*(x(il+start)-x(il+start-1));
    disp(['X1: ' num2str(X1*2)]);
    
    plot(X1,0,'ko','MarkerSize',10)
end

utop = u(:,end);
il = find(max(utop,0)==0,1);
if (~isempty(il))
    X2(j) = x(il) + (-utop(il-1)/(utop(il)-utop(il-1)))*(x(il+1)-x(il));
    disp(['X2: ' num2str(X2*2)]);
    plot(X2,0,'ko','MarkerSize',10)
end
il = find(min(utop(il:end),0)==0,1) + il - 1;
if (~isempty(il))
    X3(j) = x(il) + (-utop(il-1)/(utop(il)-utop(il-1)))*(x(il+1)-x(il));
    disp(['X3: ' num2str(X3*2)]);
    plot(X3,0,'ko','MarkerSize',10)
end

%% velocity
velocity

%% profiles at x=7 and x=15
% figure
% hold on
% i7 = find(x==7);
% plot(u(i7-1,:),yp,'kx-');
% i15=find(x==15);
% plot(u(i15-1,:),yp,'bx-');
% 
% y_Gartling = [0.5 0.45 0.4 0.35 0.3 0.25 0.20 0.15 0.10 0.05 0.0 -0.05 -0.10 -0.15 -0.2 -0.25 -0.3 -0.35 -0.4 -0.45 -0.5];
% u_Gartling7 = [0.0 -0.038 -0.049 -0.032 0.015 0.092 0.204 0.349 0.522 0.709 0.885 1.024 1.105 1.118 1.062 0.948 0.792 0.613 0.428 0.232 0.000];
% u_Gartling15 = [0.0 0.101 0.202 0.304 0.408 0.512 0.613 0.704 0.779 0.831 0.853 0.844 0.804 0.737 0.649 0.547 0.438 0.328 0.218 0.109 0.00];
% 
% % u_ex15 = (a/8)*(yp+L_y/2).*(L_y/2-yp);
% 
% plot(u_Gartling7,y_Gartling,'ro-');
% plot(u_Gartling15,y_Gartling,'go-');