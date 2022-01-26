clear all
% profile on
% N = 40000;
% s = sparse(N,1);
% % s([1:40]*100,1) = 1;
% s([1:2]*100,1) = 1;
% s.*s;
% spdiags(s,0,N,N)*s;
% diag(s)*s;
% profile off
% profile viewer

%% 
% a = (1:5)';
% 
% func(1).*a;
% disp('now')
% diag(func(1))*a;
% 
% function res = func(x)
% disp('func')
% res = x;
% end

%% 
profile on
N = 40000;
a = (1:N)';
A = spdiags(a,0,N,N);
a.*A;
A.*a;
spdiags(a,0,N,N)*A;
profile off
profile viewer

%%
% profile on
% N = 40000;
% a = 4;
% A = 1:N^2;
% a.*A;
% a*A;
% profile off
% profile viewer




