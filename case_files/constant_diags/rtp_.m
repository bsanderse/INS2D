figure(2)

[up,vp,qp] = get_velocity(V,t,options);

subplot(2,2,1)
contour(xp,yp,up')
axis square
colorbar

title("u velocity") 


subplot(2,2,2)
contour(xp,yp,vp')
axis square
colorbar

title("v velocity")

subplot(2,2,3)
plot(xp,up)
axis square

title("u velocity slices")

subplot(2,2,4)

title("time = " + num2str(t))
