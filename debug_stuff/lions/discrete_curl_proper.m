clear all;
% periodic 2x2 grid

curl = -[-1 0 1 0 1 -1 0 0; 0 -1 0 1 -1 1 0 0; 
        1 0 -1 0 0 0 1 -1; 0 1 0 -1 0 0 -1 1];
    
Dx = [-1 1 0 0; 1 -1 0 0; 0 0 -1 1; 0 0 1 -1];

Dy = [-1 0 1 0; 0 -1 0 1; 1 0 -1 0; 0 1 0 -1];

% cross1 = [ zeros(4) [zeros(2) eye(2); eye(2) zeros(2)]; -eye(4) zeros(4)];
cross1 = [ zeros(4) eye(4); -eye(4) zeros(4)];

helper = [eye(4) eye(4)];

Dxy = blkdiag(Dx,Dy);

% helper*Dxy*cross1 - curl

%% verified: helper*Dxy*cross1 = curl
if norm(helper*Dxy*cross1 - curl) ~= 0
    error('basics fail')
end
%%

psi = [2;4;8;16];
F = [3 9 27 81 1 5 25 125]';

psi_d = diag(psi);
psi_2 = [psi;psi];
psi_d2 = diag([psi; psi]);

LHS1 = curl*psi_d2*F

LHS2 = helper*Dxy*cross1*psi_d2*F

RHS1 = diag(psi)*curl*F

diff = curl*psi_d2*F-diag(psi)*curl*F
% 
% helper*diag(Dxy*(psi_2))
% 
RHS2 = helper*diag(Dxy*(psi_2))*cross1*F

helper2 = [eye(4) [0 0 0 1; 0 0 1 0; 0 1 0 0; 1 0 0 0]];
helper3 = - helper2([2 1 4 3],:)

RHS2 = helper3*diag(Dxy*(psi_2))*cross1*F

LHS1 - (RHS1+RHS2)

% interesting: helper*Dxy*psi_d2*cross1*F == curl*psi_d2*F 
%           == helper*Dxy*cross1*psi_d2*F
if norm(helper*Dxy*psi_d2*cross1*F - curl*psi_d2*F) > 0
    warning('mistake!')
end



