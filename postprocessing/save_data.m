%% save only important data
% filename = ['Ek1_N',num2str(Nx),'_t',num2str(t_end)];

if (strcmp(visc,'LES'))
    
save(filename,'uh','vh','p','up','vp','qp','omega_temp','x','y','xin','yin','xp','yp',...
     'Npx','Npy','Omu','Omv','dt','t_end','Re',...
     'convu','convv','freq','Ek1','Ek2','nu_t')
 
else
  
save(filename,'uh','vh','p','up','vp','qp','omega_temp','x','y','xin','yin','xp','yp',...
     'Npx','Npy','Omu','Omv','dt','t_end','freq','Ek1','Ek2','Re',...
     'convu','convv')
 
end
    