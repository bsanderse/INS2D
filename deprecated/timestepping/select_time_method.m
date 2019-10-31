if (method==2)
  time_AB_CN;
elseif (method==5)
  time_oneleg;
elseif (method==20)
  time_ERK;
elseif (method==21)
  time_IRK;
else
  error('time integration method unknown');
end