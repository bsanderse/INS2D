nzeros = filelen-length(num2str(n));
n_new  = num2str(10^nzeros);
file_restart = [path_results '/restart_'  n_new(2:end) num2str(n) '.mat'];
save (file_restart,'uh','vh','p','n','t','dt');