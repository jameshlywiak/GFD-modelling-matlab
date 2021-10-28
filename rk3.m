function u = rk3(u,c,dx,order,bcflag,dt)
%
% time steps the variable u through a single step of size dt
% usage is
% u = rk3(u,c,dx,order,bcflag,dt,eqflag)
% where u is the variable at the present time level
%       u is the variable at the next time level
%       dt is the time step
%       c is the constant advection speed (case1) OR coefficient array (case 2,3,4)
%       bcflag is the boundary condition
%       order is the order of the differencing scheme
%

    r = rhsadv(u,c,dx,order,bcflag);
    ut = u + dt*r;

    r = rhsadv(ut,c,dx,order,bcflag);
    ut = 0.75*u + 0.25*(ut + dt*r);

    r = rhsadv(ut,c,dx,order,bcflag);
    u = (u + 2.0*(ut + dt*r))/3.0;

end