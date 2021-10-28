function [ux0,ux,xe0,xe,Err,dx] = FDA5th(N)
% Npts calculates the exact derivative of function u0 and a 5th order 
%   finite difference approximation for a desired grid size N+1
%   *Code specific to my CFD hw*

% Exact solution

xmin0=-1; xmax0=1;
[xe0,dx0]=FDGrid(xmin0,xmax0,N);

u0=zeros(1,length(xe0));
ux0=zeros(1,length(xe0));
alpha=4;

for i=1:length(xe0)
    
    u0(i)=tanh(alpha.*(cos(pi.*xe0(i))+.5))+exp(-(cos(pi.*xe0(i)).^2));
    ux0(i)=pi.*sin(pi.*xe0(i)).*(alpha.*((tanh(alpha.*(cos(pi.*xe0(i))+.5))).^2-1)+...
        2*cos(pi.*xe0(i)).*exp(-(cos(pi.*xe0(i)).^2)));
    
end

% 5th order Scheme

% Define the domain and the computational grid;

xmin=-1; xmax=1;
[xe,dx]=FDGrid(xmin,xmax,N);

% Define initial condition
u=zeros(1,length(xe));                    % initial function
alpha=4;

for i=1:length(xe)
    
    u(i)=tanh(alpha*(cos(pi*xe(i))+.5))+exp(-(cos(pi*xe(i)).^2));

end

% Compute ux
ux=zeros(1,length(xe));
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

% 

for i=(abs(l)+1):(length(u)-r)

    ux(i)=(1/(60*dx)).*(a1*u(i+l)+a2*u(i+l+1)+a3*u(i+l+2)+a4*u(i+l+3)...
        +a5*u(i+l+4)+a6*u(i+l+5));
     
end

% Error calculation

Err=(ux0(:)-ux(:));

end

