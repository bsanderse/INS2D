% test cg with general (sparse) positive definite matrix
clear all
clc
format long
format compact

addpath('~/Software/petsc-3.1-p4/bin/matlab');

n = 10;
N = 500000;

acc   = 1e-14;
maxit = 1000;

diag1 = ones(N,1);
diags =[-2;-1;0;1;2];
D = spdiags([diag1 diag1 diag1 diag1 diag1],diags,N,N);

A = sprand(D) + speye(N)*(N-1);
% make symmetric
C = A+A';

% for cg:
[E,d] = spdiags(C);
nd = (length(d)+1)/2;
dia = d(nd:end);
E   = E(:,nd:-1:1);

time_direct = 0;
time_cg     = 0;
time_petsc  = 0;

% open socket only once
% launch('./petsc_poisson_v2 ',1,' ');
PS = PetscOpenSocket;
PetscBinaryWrite(PS,C);

for i=1:n
    
    p = rand(N,1);

    %direct solver
    tic
    x = C\p;
    time_direct = time_direct + toc;


    tic
    [x3,iter,norm1,norm2]=cg(double(E),int64(dia),int64(nd),double(p),double(acc),int64(N),double(p),int64(maxit));
    time_cg = time_cg + toc;
    
%     norm(x3-x,'inf')/norm(x3,'inf')

    tic  
    PetscBinaryWrite(PS,p');

    x4 = PetscBinaryRead(PS);
    iter = PetscBinaryRead(PS);
    norm = PetscBinaryRead(PS);

%     norm(x4'-x,'inf')/norm(x4','inf')
    time_petsc = time_petsc + toc;


end
    close(PS);

time_direct
time_cg
time_petsc