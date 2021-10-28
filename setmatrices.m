function [fp,fu,fv,fz] = setmatrices(M,N)
%
% Define the arrays for the tests
% Usage is  [fp,fu,fv,fz] = TestCgridSetMatrices(M,N)
% M-N is the number of cells in the x-y directions
%

  fp = zeros([M   N  ]);
  fu = zeros([M+1 N  ]);
  fv = zeros([M   N+1]);
  fz = zeros([M+1 N+1]);
  for j = 1:N
    fp(1:M,j) = (2*j-1)*(2*M+1) + [1:2:2*M-1]';
  end
  for j = 1:N+1
    fv(1:M,j) = 2*(j-1)*(2*M+1) + [1:2:2*M-1]';
  end
  for j = 1:N
    fu(1:M+1,j) = (2*j-1)*(2*M+1) + [0:2:2*M]';
  end
  for j = 1:N+1
    fz(1:M+1,j) = 2*(j-1)*(2*M+1) + [0:2:2*M]';
  end

 %disp '================Original Arrays---============'
 %prettyprint4(fz,fv,fu,fp,'Original Input Matrices (fp,fu,fv,fz) Laid out on a C-grid');
