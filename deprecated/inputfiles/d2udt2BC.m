function dudt = d2udt2BC(x,y,t,Re)

% Taylor:

    dudt = (2*(pi^2)/Re)^2*sin(pi*x)*cos(pi*y)*exp(-2*(pi^2)*t/Re);

end