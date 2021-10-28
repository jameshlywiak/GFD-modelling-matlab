  M = 5;
  N = 4;
  am = [0.5  0.5];      % averaging coefficients
  ad = [1.0 -1.0];      % gradient coefficients

  [fp,fu,fv,fz]=setmatrices(M,N); % define functions on C-grid points
% size(fp)
% size(fu)
% size(fv)
% size(fz)

  gu = -9.99*ones(size(fu));
  gz = -9.99*ones(size(fz));
  gv = -9.99*ones(size(fv));
  gp = -9.99*ones(size(fp));

  disp('============x-operations with xop2============');
  gu = xop2_2d(gu, fp, 0.0, am); % p to u x-average
  prettyprint4(fz,fv,gu,fp,'x-mean operation from p to u-points')
  disp('Press enter for more test');
  pause

  gp = xop2_2d(gp, fu, 0.0, am); % u to p x-average
  prettyprint4(fz,fv,fu,gp,'x-mean operation from u to p-points')
  disp('Press enter for more test');
  pause

  gz = xop2_2d(gz, fv, 0.0, am); % v to z x-average
  prettyprint4(gz,fv,fu,fp,'x-mean operation from v to z-points')
  disp('Press enter for more test');
  pause

  gv = xop2_2d(gv, fz, 0.0, am); % z to v x-average
  prettyprint4(fz,gv,fu,fp,'x-mean operation from z to v-points')
  disp('Press enter for more test');
  pause

  gu = -9.99*ones(size(fu));
  gz = -9.99*ones(size(fz));
  gv = -9.99*ones(size(fv));
  gp = -9.99*ones(size(fp));

  disp('============y-operations with yop2============');
  gv = yop2_2d(gv, fp, 0.0, am); % p to v y-average
  prettyprint4(fz,gv,fu,fp,'y-mean operation from from p to v-points')
  disp('Press enter for more test');
  pause

  gp = yop2_2d(gp, fv, 0.0, am); % v to p y-average
  prettyprint4(fz,fv,fu,gp,'y-mean operation from from v to p-points')
  disp('Press enter for more test');
  pause

  gz = yop2_2d(gz, fu, 0.0, am); % u to z y-average
  prettyprint4(gz,fv,fu,fp,'y-mean operation from from u to z-points')
  disp('Press enter for more test');
  pause

  gu = yop2_2d(gu, fz, 0.0, am); % z to u y-average
  prettyprint4(fz,fv,gu,fp,'y-mean operation from from z to u-points')
