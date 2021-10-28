function ux = FDA1(u,dx)
%Calculates 1st order FDA of function u

ux=zeros(1,length(u));

% periodic bc

ux(1)=(u(1)-u(end-1))/dx;

% compute interior

for i=2:length(ux)

    ux(i)=(u(i)-u(i-1))/dx;
     
end

end

