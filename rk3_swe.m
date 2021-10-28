function [u,v,h] = rk3_swe(u,v,h,dx,dt,H,g,f,M,N)
% RK3 timestep scheme for the non-linear SWE
% case 1,2,3 -> operation in x-direction, y-direction, and both

if size(H,1)==1
    U = H*u; V = H*v; Iota = zeros(M,N);
    Iota = g*h+((u(2:M+1,:).^2+u(1:M,:).^2)+...
        v(:,2:N+1).^2+v(:,1:N).^2)/4;
else
    Hx = (H(2:end,:)+H(1:end-1,:))/2; Hy = (H(:,2:end)+H(:,1:end-1))/2;
    U = Hx(:,1:end-2).*u; V = Hy(1:end-2,:).*v; 
    Iota = g*h+((u(2:end,:).^2+u(1:end-1,:).^2)+...
        v(:,2:end).^2+v(:,1:end-1).^2)/4;
end

ru = -xop2_2d(u,Iota,0,[1/dx -1/dx]);
ru = ru + nonlinear_calc(u,v,V,H,f,dx,M,N,1);
rv = -yop2_2d(v,Iota,0,[1/dx -1/dx]); 
rv = rv - nonlinear_calc(u,v,U,H,f,dx,M,N,2);
rh = -xop2_2d(h,U,0,[1/dx -1/dx]) - yop2_2d(h,V,0,[1/dx -1/dx]);
ut = u + dt*ru;
vt = v + dt*rv;
ht = h + dt*rh;

if size(H,1)==1
    U = H*ut; V = H*vt; Iota = zeros(M,N);    
    Iota = g*ht+((ut(2:M+1,:).^2+ut(1:M,:).^2)+...
        vt(:,2:N+1).^2+vt(:,1:N).^2)/4;
else
    U = Hx(:,1:end-2).*ut; V = Hy(1:end-2,:).*vt; 
    Iota = g*ht+((ut(2:end,:).^2+ut(1:end-1,:).^2)+...
        vt(:,2:end).^2+vt(:,1:end-1).^2)/4;
end

ru = -xop2_2d(ut,Iota,0,[1/dx -1/dx]);
ru = ru + nonlinear_calc(ut,vt,V,H,f,dx,M,N,1);
rv = -yop2_2d(vt,Iota,0,[1/dx -1/dx]); 
rv = rv - nonlinear_calc(ut,vt,U,H,f,dx,M,N,2);
rh = -xop2_2d(ht,U,0,[1/dx -1/dx]) - yop2_2d(ht,V,0,[1/dx -1/dx]);
ut = 0.75*u + 0.25*(ut + dt*ru);
vt = 0.75*v + 0.25*(vt + dt*rv);
ht = 0.75*h + 0.25*(ht + dt*rh);

if size(H,1)==1
    U = H*ut; V = H*vt; Iota = zeros(M,N);
    Iota = g*ht+((ut(2:M+1,:).^2+ut(1:M,:).^2)+...
        vt(:,2:N+1).^2+vt(:,1:N).^2)/4;
else
    U = Hx(:,1:end-2).*ut; V = Hy(1:end-2,:).*vt;
    Iota = g*ht+((ut(2:end,:).^2+ut(1:end-1,:).^2)+...
        vt(:,2:end).^2+vt(:,1:end-1).^2)/4;
end

ru = -xop2_2d(ut,Iota,0,[1/dx -1/dx]);
ru = ru + nonlinear_calc(ut,vt,V,H,f,dx,M,N,1);
rv = -yop2_2d(vt,Iota,0,[1/dx -1/dx]); 
rv = rv - nonlinear_calc(ut,vt,U,H,f,dx,M,N,2);
rh = -xop2_2d(ht,U,0,[1/dx -1/dx]) - yop2_2d(ht,V,0,[1/dx -1/dx]);
u = (u + 2.0*(ut + dt*ru))/3.0;
v = (v + 2.0*(vt + dt*rv))/3.0;
h = (h + 2.0*(ht + dt*rh))/3.0;

end

