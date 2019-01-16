% actuator disk
Ct = 0.7;   % thrust coefficient
D  = 1;     % diameter
 
m_turb = 5;
n_turb = 1;
d_turb = 4;
th_turb = 90;
tot_turb = m_turb*n_turb;

% assuming uniform grid
nd = D/deltay;

x0 = 5;
y0 = 5 + m_turb*d_turb*cosd(th_turb);

Fx = zeros(Nu,1);

for j_turb=1:n_turb
    for i_turb=1:m_turb
        s_turb = i_turb + (j_turb-1)*m_turb;
    
 x_turb(s_turb) = x0 + (i_turb-1)*d_turb*sind(th_turb) + (j_turb-1)*d_turb*cosd(th_turb);
 y_turb(s_turb) = y0 - (i_turb-1)*d_turb*cosd(th_turb) + (j_turb-1)*d_turb*sind(th_turb);
    
        % index of closest u-volume
        [val,xpos] = min(abs(x_turb(s_turb)-xin));
        [val,ypos] = min(abs(y_turb(s_turb)-yp));
        
        yrange     = ypos - floor(nd/2)+1:ypos + floor(nd/2);
    
        y1Dx       = zeros(Nux_in,1);
        y1Dx(xpos) = 1;
        y1Dy       = zeros(Nuy_in,1);
        y1Dy(yrange) = hy(yrange);                          % uniform loading
        % y1Dy(yrange)= const*force_distr.*hy(yrange);      % non-uniform loading
        Fx          = Fx - 0.5*Ct*kron(y1Dy,y1Dx);

        
    end
end

% plot(x_turb,y_turb,'rx');
