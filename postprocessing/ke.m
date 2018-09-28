% k and e



figure(5)
contour(xp,yp,reshape(kth,Npx,Npy)',20);
axis equal
axis([x1 x2 y1 y2]);
xlabel('x');
ylabel('y');
 
figure(6)
contour(xp,yp,reshape(eh,Npx,Npy)',20);
axis equal
axis([x1 x2 y1 y2]);
xlabel('x');
ylabel('y');
