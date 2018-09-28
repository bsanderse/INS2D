function dvdt = d2vdt2BC(x,y,t,Re)

     dvdt = (-2*(pi^2)/Re)^2*cos(pi*x)*sin(pi*y)*exp(-2*(pi^2)*t/Re);

end