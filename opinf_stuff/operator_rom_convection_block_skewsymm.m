function [conv_quad] = operator_rom_convection_block_skewsymm(options)
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

    K = blkdiag([Cux Cuy], [Cvx Cvy]);
    A = blkdiag([Au_ux; Au_uy], [Av_vx; Av_vy]);
    I = [blkdiag(Iu_ux, Iv_uy); blkdiag(Iu_vx, Iv_vy)];

    NF = size(A,1);

    % very inefficient
    % P = [];
    % for i=1:NF
    %     P = [P sparse(i,i,1,NF,NF)];
    % end

    % more efficient
    x = (1:NF)';
    y = x + NF*(x-1);
    P = sparse(x,y,1,NF,NF^2);

    C = K*P*kron(I,A);

    % even more efficient: don't construct P explicitly

    % even more efficient: replace P*kron(I,A) by Hadamard product 
    % -> ! I.*A does not work !


    %% switch to ROM operator
    Phi = options.rom.B;

    % conv_quad = Phi'*C*kron(Phi,Phi); % causes storage problems

    r = size(Phi,2);

    conv_quad = zeros(r,r,r);
    for i = 1:r
        conv_quad(:,:,i) = Phi'*C*kron(Phi(:,i),Phi);
    end

    conv_quad = conv_quad(:,:);