function ru = nonlinear_calc(u,v,U,H,f,dx,M,N,uvflag)
%This code is used to solve the RHS of the nonlinear SWE eq

q = zeros(M+1,N+1);
q(1,1:N-1) = (u(1,2:N)-u(1,1:N-1)+f(1,1:N-1)*dx)/(dx*H);        % BC
q(1:M-1,1) = (v(2:M,1)-v(1:M-1,1)+f(1:M-1,1)*dx)/(dx*H);        % BC

q(2:M,2:N) = (v(2:M,2:N)-v(1:M-1,2:N))/(H*dx)...
    - (u(2:M,2:N)-u(2:M,1:N-1))/(H*dx) + f(2:M,2:N)/H;


switch uvflag
    case 1    % solving for u using V
        
        Vavgx = zeros(M+1,N+1);
        qV = zeros(M+1,N+1);
        ru = zeros(M+1,N);
        
        Vavgx(2:M,:) = (U(2:M,1:N+1)+U(1:M-1,1:N+1))/2;
        qV = q.*Vavgx;    
        ru = (qV(:,2:N+1)+qV(:,1:N))/2;
    case 2    % solving for v using U
        
        Uavgy = zeros(M+1,N+1);
        qU = zeros(M+1,N+1);
        ru = zeros(M,N+1);
        
        Uavgy(:,2:N) = (U(1:M+1,2:N)+U(1:M+1,1:N-1))/2;
        qU = q.*Uavgy;
        ru = (qU(2:M+1,:)+qU(1:M,:))/2;
end

end

