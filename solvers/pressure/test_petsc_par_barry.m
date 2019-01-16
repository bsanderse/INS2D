% test cg with general (sparse) positive definite matrix
clear all
clc
format long
format compact

addpath('/net/shareware/src/petsc-3.1-p4/bin/matlab/');

N = 100000;

% create spd matrix
diag1 = ones(N,1);
diags =[-2;-1;0;1;2];
D = spdiags([diag1 diag1 diag1 diag1 diag1],diags,N,N);
A = sprand(D) + speye(N)*(N-1);

% make symmetric
C = A+A';

% right-hand side
p = rand(N,1);
time_petsc  = 0;

% open socket
PS = PetscOpenSocket(5600);

PetscBinaryWrite(PS,p');

x    = PetscBinaryRead(PS);

x(1)
p(1)


close(PS);
