
%% divergence

% Mx
ybc    = kron(uLe_i,Mx_BC.ybc1) + kron(uRi_i,Mx_BC.ybc2);
yMx    = Mx_BC.Bbc*ybc;
if (order4==1)
   ybc3 = kron(uLe_i,Mx_BC3.ybc1) + kron(uRi_i,Mx_BC3.ybc2);
   yMx3 = Mx_BC3.Bbc*ybc3;
   yMx  = alfa*yMx - yMx3;
end

% My 
ybc    = kron(My_BC.ybc1,vLo_i) + kron(My_BC.ybc2,vUp_i);
yMy    = My_BC.Bbc*ybc;
if (order4==1)
    ybc3 = kron(My_BC3.ybc1,vLo_i) + kron(My_BC3.ybc2,vUp_i);
    yMy3 = My_BC3.Bbc*ybc3;
    yMy  = alfa*yMy - yMy3;
end

yM     = yMx + yMy;

% derivative of divergence
ybc    = kron(dudtLe_i,Mx_BC.ybc1) + kron(dudtRi_i,Mx_BC.ybc2);
ydMx   = Mx_BC.Bbc*ybc;
if (order4==1)
    ybc3  = kron(dudtLe_i,Mx_BC3.ybc1) + kron(dudtRi_i,Mx_BC3.ybc2);
    ydMx3 = Mx_BC3.Bbc*ybc3;
    ydMx  = alfa*ydMx - ydMx3;
end

% My 
ybc    = kron(My_BC.ybc1,dvdtLo_i) + kron(My_BC.ybc2,dvdtUp_i);
ydMy   = My_BC.Bbc*ybc;
if (order4==1)
    ybc3  = kron(My_BC3.ybc1,dvdtLo_i) + kron(My_BC3.ybc2,dvdtUp_i);
    ydMy3 = My_BC3.Bbc*ybc3;
    ydMy  = alfa*ydMy - ydMy3;
end

ydM    = ydMx + ydMy;

% if (ibm==1)
%     ydM = [ydM; zeros(n_ibm,1)];
%     yM = [yM; zeros(n_ibm,1)];
% end