function [C_hom,C_hom_inhom,C_hom_bc,C_inhom,C_inhom_bc,C_bc] = ...
    operator_rom_convection_unsteadyBC(P,options)

Cux = options.discretization.Cux;
Cuy = options.discretization.Cuy;
Cvx = options.discretization.Cvx;
Cvy = options.discretization.Cvy;

Au_ux = options.discretization.Au_ux;
Au_uy = options.discretization.Au_uy;
Av_vx = options.discretization.Av_vx;
Av_vy = options.discretization.Av_vy;

Iu_ux = options.discretization.Iu_ux;
Iv_uy = options.discretization.Iv_uy;
Iu_vx = options.discretization.Iu_vx;
Iv_vy = options.discretization.Iv_vy;

A = blkdiag([Au_ux; Au_uy],[Av_vx; Av_vy]); %not efficient but intuitive
I = [blkdiag(Iu_ux,Iv_uy);blkdiag(Iu_vx,Iv_vy)]; %not efficient but intuitive
K = blkdiag([Cux Cuy],[Cvx Cvy]); %not efficient but intuitive

B = options.rom.B;
M1  = size(P,1);  
M   = options.rom.M;
Vbc0 = options.rom.Vbc0;
M_inhom = size(Vbc0,2);

yAu_ux = options.bc_options1.discretization.yAu_ux;
yAu_uy = options.bc_options1.discretization.yAu_uy;
yAv_vx = options.bc_options1.discretization.yAv_vx;
yAv_vy = options.bc_options1.discretization.yAv_vy;

yIu_ux = options.bc_options1.discretization.yIu_ux;
yIv_uy = options.bc_options1.discretization.yIv_uy;
yIu_vx = options.bc_options1.discretization.yIu_vx;
yIv_vy = options.bc_options1.discretization.yIv_vy;

y_A1 = [yAu_ux;yAu_uy;yAv_vx;yAv_vy];
y_I1 = [yIu_ux;yIv_uy;yIu_vx;yIv_vy];

yAu_ux = options.bc_options2.discretization.yAu_ux;
yAu_uy = options.bc_options2.discretization.yAu_uy;
yAv_vx = options.bc_options2.discretization.yAv_vx;
yAv_vy = options.bc_options2.discretization.yAv_vy;

yIu_ux = options.bc_options2.discretization.yIu_ux;
yIv_uy = options.bc_options2.discretization.yIv_uy;
yIu_vx = options.bc_options2.discretization.yIu_vx;
yIv_vy = options.bc_options2.discretization.yIv_vy;

y_A2 = [yAu_ux;yAu_uy;yAv_vx;yAv_vy];
y_I2 = [yIu_ux;yIv_uy;yIu_vx;yIv_vy];

y_A = [y_A1 y_A2];
y_I = [y_I1 y_I2];
M_bc = size(y_A,2); 

C_hom = zeros(M1,M*M);
for i=1:M
    for j=1:M
        Conv = K*((I*B(:,i)).*(A*B(:,j)));
        % convert third order tensor to a matrix via the following indexing:
        k = sub2ind([M,M],j,i); % first loop over j,  then over i
        C_hom(:,k) = P*Conv;        
    end
end

% observe that this implementation is different to the theoretical
% description in the Master thesis
C_hom_inhom = zeros(M1,M*M_inhom);
for i=1:M
    for j=1:M_inhom
        Conv = K*((I*B(:,i)).*(A*Vbc0(:,j))) + K*((I*Vbc0(:,j)).*(A*B(:,i)));
        % convert third order tensor to a matrix via the following indexing:
        k = sub2ind([M_inhom,M],j,i); % first loop over j,  then over i
        C_hom_inhom(:,k) = P*Conv;        
    end
end

C_hom_bc = zeros(M1,M*M_bc);
for i=1:M
    for j=1:M_bc
        Conv = K*((I*B(:,i)).*y_A(:,j)) + K*((y_I(:,j)).*(A*B(:,i)));
        % convert third order tensor to a matrix via the following indexing:
        k = sub2ind([M_bc,M],j,i); % first loop over j,  then over i
        C_hom_bc(:,k) = P*Conv;        
    end
end

C_inhom = zeros(M1,M_inhom*M_inhom);
for i=1:M_inhom
    for j=1:M_inhom
        Conv = K*((I*Vbc0(:,i)).*(A*Vbc0(:,j)));
        % convert third order tensor to a matrix via the following indexing:
        k = sub2ind([M_inhom,M_inhom],j,i); % first loop over j,  then over i
        C_inhom(:,k) = P*Conv;        
    end
end

C_inhom_bc = zeros(M1,M_inhom*M_bc);
for i=1:M_inhom
    for j=1:M_bc
        Conv = K*((I*Vbc0(:,i)).*y_A(:,j)) + K*((y_I(:,j)).*(A*Vbc0(:,i)));
        % convert third order tensor to a matrix via the following indexing:
        k = sub2ind([M_bc,M_inhom],j,i); % first loop over j,  then over i
        C_inhom_bc(:,k) = P*Conv;        
    end
end

C_bc = zeros(M1,M_bc*M_bc);
for i=1:M_bc
    for j=1:M_bc
        Conv = K*((y_I(:,i)).*(y_A(:,j)));
        % convert third order tensor to a matrix via the following indexing:
        k = sub2ind([M_bc,M_bc],j,i); % first loop over j,  then over i
        C_bc(:,k) = P*Conv;        
    end
end
