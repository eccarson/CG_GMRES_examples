% CG implementation based on 3-term recurrences (STCG)
% Uses simulated exact arithmetic
% Note: Requires Advanpix toolbox
%
% Input:
%   A: coefficient matrix
%   b: right-hand side in Ax=b
%   x0: initial approximate solution
%   nmax: maximum number of iterations
%   nd: number of digits to use to simulate exact arithmetic
%
% Output:
%   x: approximate solution
%   r: residual vector
%   Aerr: A-norm of the error

function [ x,r,Aerr] = STCG_mp( A,b,x0,nmax,nd)

mp.Digits(nd);

b = mp(b, nd);
A = mp(A, nd);

xtrue = A\b;
x(:,1) = mp(x0,nd);
r0=mp(b-A*x(:,1),nd);
r(:,1) = r0;

err(:,1) = mp(xtrue - x(:,1),nd);
Aerr(1) = mp(sqrt(mp(err(:,1)'*A*err(:,1)))/sqrt(xtrue'*A*xtrue));

for i=1:nmax
    if i==1
        eprev = mp(0);
        xprev = mp(x0);
        rprev = mp(r0);
    else
        eprev = e(i-1);
        xprev = x(:,i-1);
        rprev = r(:,i-1);
    end
    q(i) = mp((r(:,i)'*A*r(:,i))/(r(:,i)'*r(:,i))-eprev,nd);
    x(:,i+1) = mp(x(:,i) + (1/q(i))*(r(:,i) + eprev*(x(:,i)-xprev)),nd);
    r(:,i+1) = mp(r(:,i) + (1/q(i))*(-A*r(:,i)+eprev*(r(:,i)-rprev) ),nd);
    e(i) = mp(q(i)* ((r(:,i+1)'*r(:,i+1))/(r(:,i)'*r(:,i))),nd);

    err(:,i+1) = mp(xtrue - x(:,i+1),nd);
    Aerr(i+1) = mp(sqrt(mp(err(:,i+1)'*A*err(:,i+1)))/sqrt(xtrue'*A*xtrue));
    
end


end

