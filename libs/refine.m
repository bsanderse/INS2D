function [xnew] = refine(x)

xnew = zeros(1,length(x)*2-1);

j=1;
for i=1:length(x)-1
    
    xnew(j)   = x(i);
    xnew(j+1) = (x(i)+x(i+1))/2;
   
    j = j+2;
end
xnew(end) = x(end);