function r = FEuler(rhs,dx,xyflag,m,n)
% Differencing scheme for CFD hw 4, problem 3
% pflag = 0,1 for x,y

r = zeros(m,n);

%disp(size(r))

switch xyflag
    case 0 % operation in x-dir
        r(1,:) = rhs(1,:); r(m,:) = rhs(m,:);
        r(2:m,:) = (rhs(2:m,:)-rhs(1:m-1,:))/dx;
    case 1 % operation in y-dir
        r(:,1) = rhs(:,1); r(:,n) = rhs(:,n);
        r(:,2:n) = (rhs(:,2:n)-rhs(:,1:n-1))/dx;     
end

r = -r;

end

