function [R,options] = ROM_initialization(V,t,options)
%% initialize reduced order solution
% we expand the part of the solution vector that is div-free in terms of
% B*R
% V = B*R + Vbc

if options.rom.bc_recon ~= 5
% get the coefficients of the ROM
R = getROM_velocity(V,t,options);

% for projected-divergence-free ROM, enforce projected divergence-freeness
% if options.rom.bc_recon == 5
else
    % construct ROM divergence operator
        Bp = options.rom.Bp;
        B = options.rom.B;
        hatM = Bp'*options.discretization.M*B;
        cond(options.discretization.M*B)
        options.rom.hatM = hatM;
        hatL = -hatM*hatM';
        options.rom.hatL = hatL;
%         condition3(j) = cond(hatL)
        condition3 = cond(hatL)
        
%         hatM2 = options.rom.Bp2'*options.discretization.M*B;
%         hatL2 = -hatM2*hatM2';
%         condition4(j) = cond(hatL2)
%         
%         hatM3 = options.rom.Bp3'*options.discretization.M*B;
%         hatL3 = -hatM3*hatM3';
%         condition5(j) = cond(hatL3)

%     Bp = options.rom.Bp;
%     hatM = options.rom.hatM;
    hatG = -hatM';
    hatL = hatM*hatG;
    
%     yM = -options.discretization.yM; % wrong! we need BC approximation
    F_M = options.discretization.F_M;
    phi_bc = options.rom.phi_bc;
    yM = - F_M*phi_bc*get_a_bc(t,options);
%     yM = F_M*phi_bc*get_a_bc(t,options);
    hatyM = Bp'*yM;
    %     phi_bc = options.rom.phi_bc;
    %     yM = phi_bc*phi_bc'*yM;

    %     bstar = hatL\(hatM*R-Bp'*yM);
    %     Rstar = R-hatG*bstar;
    %     R = Rstar;

    % QR based method
    Om = options.grid.Om;
%     [Q_,R_] = qr(hatM');
%     Q_1 = Q_(:,1:Mp);
%     Q_2 = Q_(:,Mp+1:end);
%     R_1 = R_(1:Mp,1:Mp);
%     R = Q_1*hatyM/(R_1') + Q_2*Q_2'*B'*(Om.*V);
% 
%         %test
%     norm(hatM*R-hatyM)

    %% second try
%     M_h = options.discretization.M;
%     [Q__,R__] = qr((M_h*B)');
%     H = rank(M_h*B);
%     Q_1_ = Q__(:,1:H);
%     Q_2_ = Q__(:,H+1:end);
%     R_1_ = R__(1:H,:);
%     a_1 = (R_1_*R_1_')\(R_1_*yM);
%     a_2 = Q_2_'*B'*(Om.*V);
%     R = Q_1_*a_1 + Q_2_*a_2;
% 
%     % tests
%     norm(M_h*B*Q_1_*a_1-yM)
%     norm(M_h*B*R - yM)
%     diff_ = B*R-V;
%     norm(diff_)
%     norm(M_h*diff_)
%     norm(a_2)
%     norm(diff_'*(Om.*V))
%     diff_2 = B*Q_1_*a_1 - V;

%botch
% qr_test
    


%% botch
    
    if options.rom.bases_construction == "mthesis" || ...
            options.rom.bases_construction == "optimal"
        phi_hom = options.rom.phi_hom;
        ahom = phi_hom'*(Om.*V);
        
        % compare two methods to compute ahom -> indeed equivalent
%         ahom1 = phi_hom'*(Om.*V);

%         V = V - get_unsteadyVbc(t,options); %not sure if necessary as Vbc should be orthogonal to B
%         ahom = phi_hom'*(Om.*V);
% 
%         norm(ahom1-ahom)

        ainhom = get_a_inhom(t,options);
        
        R = [ahom(:); ainhom(:)];
    elseif options.rom.bases_construction == "qr"
            %% new try
    M_h = options.discretization.M;
    themat = (M_h*B)';
    [Q,R,P] = qr(themat,0);
    norm(themat(:,P)-Q*R)
    RR = R(:,P);
    norm(themat-Q*RR)
    H = rank(themat);
    Q1 = Q(:,1:H);
    Q2 = Q(:,H+1:end);
    R1 = RR(1:H,:);
    norm(themat-Q1*R1)
    norm(M_h*B*Q2)
    
    a1 = (R1*R1')\(R1*yM);
    norm(M_h*B*Q1*a1 - yM)
    
    a2 = Q2'*(B'*(Om.*V));
    
    R = Q1*a1 + Q2*a2;
    norm(M_h*B*R - yM)
    
%     R_ = Q'*(B'*(Om.*V));
%     norm(R_-R) %not 0, but why?


    end
    
    

    
    %%
%     [QQ,RR] = qr(full(M_h)); % expensive!
%     HH = rank(M_h);
%     QQ1 = QQ(:,1:HH);

%     Q1t = options.discretization.Q1t;
%     Q2t = options.discretization.Q2t;
%     R1 = options.discretization.R1;
% 
%     a1 = (R1*R1')\(R1*yM); % only valid if ROM basis includes all relevant inhomogeneous modes
%     a11 = Q1t'*(Om.*V); % only valid if ROM basis includes all relevant inhomogeneous modes
% 
%     norm(a1-a11)
% 
%     a0 = B'*(Om.*V);
%     norm(M_h*B*a0)
%     
%     Qt = [Q1t Q2t];
%     T = Qt'*(Om.*B);
%     a2 = Q2t'*(Om.*V);
%     at = [a1; a2];
%     a = (T'*T)\(T'*at); % does not yield exact solution to T*a = at, but approximation
%     norm(T*a-at)
%     norm(R-a)
%     R = a;

%     Vr = B*a;
%     norm(M_h*Vr-yM)
%     norm(Vr-V)
end