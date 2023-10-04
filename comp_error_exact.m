%Computes the A-norm of the error vector x, where X is the true solution
%vector, in simulated exact arithmetic 
%
% Note: Requires Advanpix library
%
function [err] = comp_error_exact(A,x,X)

% 200 digits of precision to simulate exact arithmetic
mp.Digits(200);

A = mp(A);
x = mp(x);
X = mp(X);

for i = 1:size(X,2)-1
   Xk = X(:,i)-x;
   err(i) = sqrt(Xk'*A*Xk);
end