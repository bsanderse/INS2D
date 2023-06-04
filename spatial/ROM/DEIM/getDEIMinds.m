function [ps,PTUinv] = getDEIMinds(us,DEIMdim)

n = size(us,1);

[~,p1] = max(abs(us(:,1)));

U = us(:,1);
P_ = @(pl) sparse(pl,1,1,n,1);
P = P_(p1);
ps = p1;

%% only visualization
% plot(U)
% hold on
% plot(p1,U(p1),'rx')
% xlim([800 900])
% legend('show')
%%

for l = 2:1:DEIMdim  %8% size(us,2)
    ul = us(:,l);
    c = (P'*U)\(P'*ul);
    r = ul - U*c;
    [~,pl] = max(abs(r));

%% only visualization
% plot(r)
% plot(pl,r(pl),'rx')

%%

    U = [U ul];
    P = [P P_(pl)];
    ps = [ps pl];
end

cond(P'*U)
PTUinv = inv(P'*U);

end

