function [u1,u2,p]=Navier_stokes(theta,t,u1,u2,p,m,nu,Pb,Tb,Nb,N1,elm_bry)
% unsteady ns equation(rectangle domain)
% u_t-v \Delta u + u.grad u+grad p = f, in \Omega 
% div u = 0, in \Omega 
% u = g, on \partial \Omega 

Tmax = t(end);

[u1,u2,p] = ns_pde_solver(theta,Tmax,m,nu,u1,u2,p,...
                        Pb,Tb,Nb,N1,elm_bry);
return
end







