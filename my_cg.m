function [Anorm] = my_cg(A,b,x0,x,tol,MAXIT)

%
% Implementation of CG that requires the solution as input
% and returns the relative A-norm of the error as output.
%
% Joerg Liesen, September 27, 2011
%

% Initialise
xnew = x0; 
rnew = b-A*xnew; nr0 = norm(rnew);
pnew = rnew;
ne0 = sqrt((x-x0)'*A*(x-x0));
Anorm(1) = 1;

for j=1:MAXIT
    Apnew = A*pnew;
    alpha = norm(rnew)^2 / (pnew'*Apnew);
    xnew = xnew + alpha*pnew;
    err = x-xnew;
    Anorm(j+1) = sqrt(err'*A*err)/ne0;
    rold = rnew;
    rnew = rold - alpha*Apnew;
    % convergence check; use A-norm of the error
    if (Anorm) < tol
        break
    end
    omega = (norm(rnew)/norm(rold))^2;
    pnew = rnew + omega*pnew;
end

