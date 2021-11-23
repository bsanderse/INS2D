% assuming 30 runs (6 different mode numbers with 5 repetitions each
%  Ms = [2 5 10 20 50 100];
%  Ms =[2 5 10 20 40]

taylor = @(x) sum(reshape(flip(x),[6 5]),2)/5;
Ms = taylor(M_list)
svd_end2 = taylor(svd_end);
precompute_end2 = taylor(precompute_end);
time_loop2 = taylor(time_loop);

%  Ms = sum(reshape(flip(M_list),[5 5]),2)/5;
%  svd_end2 = sum(reshape(flip(svd_end),[5 5]),2)/5;
%  time_loop2 = sum(reshape(flip(time_loop),[5 5]),2)/5;
%  precompute_end2 = sum(reshape(flip(precompute_end),[5 5]),2)/5; 
 
% svd_end2 = sum(reshape(svd_end,[5 6]),1);
% precompute_end2 = sum(reshape(precompute_end,[5 6]),1);
% time_loop2 = sum(reshape(time_loop,[5 6]),1);

%  svd_end2 = svd_end;
%  precompute_end2 = precompute_end;
%  time_loop2 = time_loop;

figure
loglog(Ms,svd_end2,'-s','displayname','offline - SVD')
hold on
loglog(Ms,precompute_end2,'-s','displayname','offline - precomputing operators')
loglog(Ms,time_loop2,'-s','displayname','online')
ylabel('CPU time [s]')
xlabel('number of modes')

% FOM_time = 78.519914999999997; % actuator ROM
FOM_time = 23.423704640000000; % actuator unsteady ROM 200x80
loglog([1 100],ones(1,2)*FOM_time,'displayname','FOM')
legend('show')