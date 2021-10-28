function [ux] = FDA3(u,dx)
% FDA5 calculates a 5th order finite difference approximation 
%   *Code specific to my CFD hw*

ux=zeros(1,length(u));
l=-2; r=1;                                 % Define stencil
a1=2; a0=3; a_1=-6; a_2=1;                 % Coefficients

% periodic bc

ux(1)=(a1.*u(2)+a0.*u(1)+a_1.*u(end-1)+a_2.*u(end-2))/(6*dx);
ux(2)=(a1.*u(3)+a0.*u(2)+a_1.*u(1)+a_2.*u(end-1))/(6*dx);
ux(end)=(a1.*u(2)+a0.*u(end)+a_1.*u(end-1)+a_2.*u(end-2))/(6*dx);

% compute interior

for i=(abs(l)+1):(length(u)-r)

    ux(i)=(a1*u(i+1)+a0*u(i)+a_1*u(i-1)+a_2*u(i-2))/(6*dx);
     
end

end

