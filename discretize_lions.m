clear all;

Lx = 1;
Ly = 1;

% Ns = [10 100 1000 10000 100000]; % last entry too large
Ns = [10 100 1000 10000];
% Ns = [4 7 10];

norms = zeros(numel(Ns),1);

for i = 1:numel(Ns)
    N = Ns(i);

    Nx = N-1;
    Ny = N-1;
    
    x = Lx/N;
    y = Ly/N;

    Mp = create_Mp(Nx,Ny,x,y);
    Mp_inv = inv(Mp);
%     norms(i) = norm(full(Mp_inv))
    norms(i) = svds(Mp_inv,1)
end


Mp = create_Mp(2,2,1,1);
full(Mp)

function Mp = create_Mp(Nx,Ny,x,y)

% basemat = @(N) [ zeros(1,N); eye(N) ] - [ eye(N); zeros(1,N) ];
% 
% mat_x = basemat(x);


xs = x*ones(Nx,1);
ys = y*ones(Ny,1);

vec = [-ys; -xs; ys; xs];
N = numel(vec);

Mp = spdiags([vec -vec], [0 1], N, N);

end