%% boundary condition on body
% Dirichlet value for velocity on body
uB = 0;      
vB = 0;
% Neumann value for pressure
pB = 0;


%% body description (should be a closed contour)
body;


%% find interface points with raytracing

% interface: n x 2 matrix, giving for each interface point the (i,j) indices
% inside: Nx x Ny matrix, giving for each (i,j) point indication if inside or
% not
disp('finding inside and interface points with raytracing...');

%u-volumes; inside excludes interface points
[interface_u, inside_incl_interface_u, inside_u] = find_insidepoints2D(xin,yp,x_k,y_k);
%v-volumes; inside excludes interface points
[interface_v, inside_incl_interface_v, inside_v] = find_insidepoints2D(xp,yin,x_k,y_k);
%p-volumes; inside includes interface points
[interface_p, inside_p, inside_excl_interface_p] = find_insidepoints2D(xp,yp,x_k,y_k);


%% find body intercepts and interface points
disp('finding body intercept and image points...');
[xbi_u,ybi_u,xip_u,yip_u,panel_bi_u] = find_intercept_points(interface_u,xin,yp,x_k,y_k);

[xbi_v,ybi_v,xip_v,yip_v,panel_bi_v] = find_intercept_points(interface_v,xp,yin,x_k,y_k);

[xbi_p,ybi_p,xip_p,yip_p,panel_bi_p] = find_intercept_points(interface_p,xp,yp,x_k,y_k);


%% find bilinear interpolation points and matrices
disp('setting up interpolation matrices...');
[bilin_mat_inv_u, bilin_points_u ] = set_interpolation_matrix(xbi_u,ybi_u,xip_u,yip_u,xin,yp,interface_u);
% 
[bilin_mat_inv_v, bilin_points_v ] = set_interpolation_matrix(xbi_v,ybi_v,xip_v,yip_v,xp,yin,interface_v);

[bilin_mat_inv_p, bilin_points_p ] = set_interpolation_matrix(xbi_p,ybi_p,xip_p,yip_p,xp,yp,interface_p,1,x_k,y_k,panel_bi_p);


%% determine inside and interface indices

% number of inside points
n_inside_p = sum(inside_p(:)); % inside points, including interface
n_inside_u = sum(inside_u(:)); % inside points, excluding interface
n_inside_v = sum(inside_v(:)); % inside points, excluding interface

% determine indices of inside points

% allocate arrays
index_p    = zeros(n_inside_p,1);
index_u    = zeros(n_inside_u,1);
index_v    = zeros(n_inside_v,1);
% counters
q_p        = 1; 
q_u        = 1; 
q_v        = 1;

% set rectangular region that encloses body to minimize overhead
imin = find(xp<min(x_k),1,'last')-2;
imax = find(xp>max(x_k),1,'first')+2; 
jmin = find(yp<min(y_k),1,'last')-2; 
jmax = find(yp>max(y_k),1,'first')+2;

if (isempty(imin))
    imin=1;
end
if (isempty(imax))
    imax=Nx;
end
if (isempty(jmin))
    jmin=1;
end
if (isempty(jmax))
    jmax=Ny;
end

for i=imin:imax
  for j=jmin:jmax
 
      if (inside_p(i,j)==1)
   
         ind          = sub2ind([Npx,Npy],i,j);
         index_p(q_p) = ind;
         q_p          = q_p + 1;
         
      end
      if (inside_u(i,j)==1)
   
         ind          = sub2ind([Nux_in,Nuy_in],i,j);
         index_u(q_u) = ind;
         q_u          = q_u + 1;
         
      end
      if (inside_v(i,j)==1)
   
         ind          = sub2ind([Nvx_in,Nvy_in],i,j);
         index_v(q_v) = ind;
         q_v          = q_v + 1;
         
      end
  end

end

n_interface_u  = length(xbi_u);
n_interface_v  = length(xbi_v);
n_interface_p  = length(xbi_p);

%% pressure related parameters
index_p_outer         = setdiff(1:Np,index_p);

interface_p_ind  = sub2ind([Npx,Npy],interface_p(:,1),interface_p(:,2));
pGC              = zeros(n_interface_p,1);

p_inner                 = zeros(n_inside_p,1);
[val interface_p_order] = sort(interface_p_ind);
[inters interface_p_inner test2]    = intersect(sort(index_p),val);

%% find betas, such that phi_IP = sum(beta_i * phi_i)
beta_u         = zeros(n_interface_u,4);
for qq=1:n_interface_u

    beta_u(qq,:) = squeeze(bilin_mat_inv_u(qq,:,:))'*[1; xip_u(qq); yip_u(qq); xip_u(qq)*yip_u(qq)];
    
end

beta_v         = zeros(n_interface_v,4);
for qq=1:n_interface_v

    beta_v(qq,:) = squeeze(bilin_mat_inv_v(qq,:,:))'*[1; xip_v(qq); yip_v(qq); xip_v(qq)*yip_v(qq)];
    
end

beta_p         = zeros(n_interface_p,4);
for qq=1:n_interface_p

    beta_p(qq,:) = squeeze(bilin_mat_inv_p(qq,:,:))'*[1; xip_p(qq); yip_p(qq); xip_p(qq)*yip_p(qq)];
    
end
% check: we need sum(beta)=1 for all interface points
% this holds only for Dirichlet BC
if ( any(abs(sum(beta_u,2)-1) > 1e-10) || ...
     any(abs(sum(beta_v,2)-1) > 1e-10))
    disp(['sum of interpolation coefficients not 1: ' ...
           num2str(max(abs(sum(beta_u,2)-1))) ' ' ...
           num2str(max(abs(sum(beta_v,2)-1)))]);
end



%% construct matrices with 1 on diagonal at interface location
interface_u_ind  = sub2ind([Nux_in,Nuy_in],interface_u(:,1),interface_u(:,2));
diag1            = zeros(Nu,1);
diag1(interface_u_ind)    = 1;
mat_u            = spdiags(diag1,0,Nu,Nu);
   
interface_v_ind  = sub2ind([Nvx_in,Nvy_in],interface_v(:,1),interface_v(:,2));
diag2            = zeros(Nv,1);
diag2(interface_v_ind)    = 1;
mat_v            = spdiags(diag2,0,Nv,Nv);

CD_interface     = [mat_u spalloc(Nu,Nv,0); spalloc(Nv,Nu,0) mat_v];  

%% construct matrices with 1 on diagonal at inside location
diag1            = zeros(Nu,1);
diag1(index_u)   = 1;
mat_u            = spdiags(diag1,0,Nu,Nu);
      
diag2            = zeros(Nv,1);
diag2(index_v)   = 1;
mat_v            = spdiags(diag2,0,Nv,Nv);

CD_inside        = [mat_u spalloc(Nu,Nv,0); spalloc(Nv,Nu,0) mat_v];  


%% for testing purposes:
% construct right hand side
% uB = 0;
% [rhs_u ] = set_interpolation_rhs(xbi_u,ybi_u,xip_u,yip_u,xin,yp,interface_u_q,u,uB);
% 
% vB = 0;
% [rhs_v ] = set_interpolation_rhs(xbi_v,ybi_v,xip_v,yip_v,xp,yin,interface_v_q,v,vB);
% 
% 
% % set ghost point values
% nq_u = length(xbi_u);
% for qq=1:nq_u
%    
%     i = interface_u_q(qq,1);
%     j = interface_u_q(qq,2);
%     coeff  = squeeze(bilin_mat_inv_u(qq,:,:))*squeeze(rhs_u(qq,:))';
%     u(i,j) = 2*uB - coeff'*[1 xip_u(qq) yip_u(qq) xip_u(qq)*yip_u(qq)]';
%     
% end
% 
% nq_v = length(xbi_v);
% for qq=1:nq_v
%     
%     i = interface_v_q(qq,1);
%     j = interface_v_q(qq,2);
%     coeff  = squeeze(bilin_mat_inv_v(qq,:,:))*squeeze(rhs_v(qq,:))';
%     v(i,j) = 2*vB - coeff'*[1 xip_v(qq) yip_v(qq) xip_v(qq)*yip_v(qq)]';     
%     
% end
% 
% pB = 0; % Neumann value, if different from zero we need to do some additional programming effort
% p=kron(ones(Npy,1)',xp);
% [rhs_p ] = set_interpolation_rhs(xbi_p,ybi_p,xip_p,yip_p,xp,yp,interface_p_q,p,pB);
% 
% nq_p = length(xbi_p);
% for qq=1:nq_p
%     
%     i = interface_p_q(qq,1);
%     j = interface_p_q(qq,2);
%     coeff  = squeeze(bilin_mat_inv_p(qq,:,:))*squeeze(rhs_p(qq,:))';
%     p(i,j) = coeff'*[1 xip_p(qq) yip_p(qq) xip_p(qq)*yip_p(qq)]';     
%     
% end
% 