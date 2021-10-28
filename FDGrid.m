function [xe,dx] = FDGrid(xmin,xmax,N)

dx = (xmax-xmin)/N;
xe = xmin:dx:xmax;

end
