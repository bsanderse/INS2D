function [x_new,y_new] = rotate_body(x,y,alfa)
% rotate a body described by (x,y) around the origin with angle alfa
% alfa in degrees, counterclockwise is positive

r  = sqrt(x.^2+y.^2);
th = atand(y./x);
th(y<=eps & x<=eps) = 0;
x_new = r.*cosd(th-alfa);
y_new = r.*sind(th-alfa);

