function [conv_bc,conv_linear,conv_quad] = operator_rom_convection_block_skewsymm2(P,options)
% precompute convective operators
% ! block-skew-symmetric and more efficient than
% operator_rom_convection_block_skewsymm
% projection with generic matrix P, size M1 x NV
% for momentum equation, P is B' or B'*Om_inv
% for Poisson equation, P is Bp'*M*Om_inv

NV = options.grid.Nu+options.grid.Nv;
Z  = zeros(NV,1);

B   = options.rom.B;
M   = options.rom.M;

if options.rom.rom_bc ~= 0
    Vbc = options.rom.Vbc;
end

% P can project to M or to Mp modes 
M1  = size(P,1);  

% elegant, but does not work because convection outputs componentwise 
% (see local function at the end)
% projected_convection = @(u,v) P*convection(u,v,0,options,0);

conv_bc = zeros(M1,1);

conv_i_bc = zeros(M1,M);
conv_bc_j = zeros(M1,M);

conv_linear = zeros(M1,M);

if options.rom.rom_bc ~= 0

    conv_bc = projected_convection(Vbc,Vbc);

    for i = 1:M
        u_i = B(:,i);

        conv_i_bc_ = projected_convection(u_i,Vbc);
        conv_bc_j_ = projected_convection(Vbc,u_i);

        conv_i_bc(:,i) = conv_i_bc_;
        conv_bc_j(:,i) = conv_bc_j_;

        conv_linear(:,i) = conv_i_bc_ + conv_bc_j_ - 2*conv_bc;
    end

end

conv_quad = zeros(M1,M,M);

for i = 1:M 
    u_i = B(:,i);
    for j = 1:M
        u_j = B(:,j);
        conv_quad(:,i,j) = projected_convection(u_i,u_j) - conv_i_bc(:,i) ...
            - conv_bc_j(:,j) + conv_bc;
    end
end

conv_quad = reshape(conv_quad,M1,M^2);

function proj_conv = projected_convection(u,v)
    [convu, convv] = convection(u,v,0,options,0);

    proj_conv = P*[convu; convv];
end

end
