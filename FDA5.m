function [ux,err] = FDA5(u,dx)
% FDA5 calculates a 5th order finite difference approximation 
%   *Code specific to my CFD hw*

% Compute ux
ux=zeros(1,length(u));
l=-3; r=2;                                 % Define stencil
a1=-2; a2=15; a3=-60; a4=20; a5=30; a6=-3; % Coefficients

% Define Periodic Boundaries

ux(1)=(1/(60*dx)).*(a1*u(end-3)+a2*u(end-2)+a3*u(end-1)+a4*u(1)...
    +a5*u(2)+a6*u(3));
ux(2)=(1/(60*dx)).*(a1*u(end-2)+a2*u(end-1)+a3*u(1)+a4*u(2)...
    +a5*u(3)+a6*u(4));
ux(3)=(1/(60*dx)).*(a1*u(end-1)+a2*u(1)+a3*u(2)+a4*u(3)...
    +a5*u(4)+a6*u(5));
ux(end-1)=(1/(60*dx)).*(a1*u(end-4)+a2*u(end-3)+a3*u(end-2)+a4*u(end-1)...
    +a5*u(1)+a6*u(2));
ux(end)=(1/(60*dx)).*(a1*u(end-3)+a2*u(end-2)+a3*u(end-1)+a4*u(1)...
    +a5*u(2)+a6*u(3));

for i=(abs(l)+1):(length(u)-r)

    ux(i)=(1/(60*dx)).*(a1*u(i+l)+a2*u(i+l+1)+a3*u(i+l+2)+a4*u(i+l+3)...
        +a5*u(i+l+4)+a6*u(i+l+5));
     
end

end

