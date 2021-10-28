function r = rhsadv(ue,c,dx,order,bcflag)
%
% Compute finite difference advection term
% usage is
% r = rhsadv(h,u,dx,order,bcflag)
% where h is the variable at the present time level
%       c is the flow velocity at cell edges
%       dx is the space size a vector of dimension d

switch order
% Evaluate derivative using one of the method described below

  case 1
    r = FDA1(ue,dx);           % Backward Euler 1st oder

  case 2
    r = cd2(ue,dx,bcflag);     % Centered 2nd order

  case 3
    r = FDA3(ue,dx);           % Upstream 3rd order

  case 4
    r = cd4(ue,dx,bcflag);     % Centered 4th order

  case 5
    r = FDA5(ue,dx);           % Upstream 5th order

end
  r =-c*r;