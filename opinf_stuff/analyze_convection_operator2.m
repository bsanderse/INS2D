
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


    P = [];
    for i=1:NF
        P = [P sparse(i,i,1,NF,NF)];
    end

    C = K*P*kron(I,A);

    %% switch to ROM operator
    Phi = options.rom.B;

    C_r = Phi'*C*kron(Phi,Phi);


    Diffu  = options.discretization.Diffu;
    Diffv  = options.discretization.Diffv;

    D = blkdiag(Diffu, Diffv);