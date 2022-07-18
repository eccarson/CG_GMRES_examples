% CG implementation based on 3-term recurrences (STCG)
% Uses double precision
% Note: Requires Advanpix toolbox (only for computing true solution)
%
% Input:
%   A: coefficient matrix
%   b: right-hand side in Ax=b
%   x0: initial approximate solution
%   nmax: maximum number of iterations
%
% Output:
%   x: approximate solution
%   r: residual vector
%   Aerr: A-norm of the error


function [ x,r,Aerr] = STCG(A,b,x0,nmax)

mp.Digits(100);
xtrue = mp(A)\mp(b);

% Ensure A,b are double
b = double(b);
A = double(A);

x(:,1) = x0;
r0=b-A*x(:,1);
r(:,1) = r0;

err(:,1) = mp(xtrue - x(:,1));
Aerr(1) = sqrt(err(:,1)'*A*err(:,1))/sqrt(xtrue'*A*xtrue);

for i=1:nmax
    if i==1
        eprev = 0;
        xprev = x0;
        rprev = r0;
    else
        eprev = e(i-1);
        xprev = x(:,i-1);
        rprev = r(:,i-1);
    end
    q(i) = (r(:,i)'*A*r(:,i))/(r(:,i)'*r(:,i))-eprev;
    x(:,i+1) = x(:,i) + (1/q(i))*(r(:,i) + eprev*(x(:,i)-xprev));
    r(:,i+1) = r(:,i) + (1/q(i))*(-A*r(:,i)+eprev*(r(:,i)-rprev) );
    e(i) = q(i)* ((r(:,i+1)'*r(:,i+1))/(r(:,i)'*r(:,i)));

    err(:,i+1) = mp(xtrue - x(:,i+1));
    Aerr(i+1) = mp(sqrt(mp(err(:,i+1)'*A*err(:,i+1)))/sqrt(xtrue'*A*xtrue));
    
end

end

