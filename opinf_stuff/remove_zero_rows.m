function less_rows = remove_zero_rows(rows)
% remove rows with all entries zeros

less_rows = rows(sum(abs(rows),2)~=0,:);
