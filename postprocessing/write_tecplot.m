% write tecplot file
Bmap  = kron(speye(Nuy_in),BMx);
up    = reshape( Bmap*(Au_ux * uh + yAu_ux), Npx, Npy);

Bmap  = kron(BMy,speye(Nvx_in));
vp    = reshape( Bmap*(Av_vy * vh + yAv_vy), Npx, Npy);


pp    = reshape(p,Npx,Npy);
Tp    = zeros(Nx,Ny);
% velocityTecplot2D('final',xp,yp,up,vp,pp,Tp);

% unix('tecplot final.plt');
zp=[0];
[xp_mesh,yp_mesh,zp_mesh] = meshgrid(xp,yp,zp);
wp = 0*up;
vtkwrite('test.vtk','structured_grid',xp_mesh,yp_mesh,zp_mesh,'vectors','vector_field',up,vp,wp,'scalars','pressure',pp);