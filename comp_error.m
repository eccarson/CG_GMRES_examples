%Computes the A-norm of the error vector x, where X is the true solution
%vector
function [err] = comp_error(A,x,X)

A = double(A);
x = double(x);
X = double(X);

for i = 1:size(X,2)-1
   Xk = X(:,i)-x;
   err(i) = sqrt(Xk'*A*Xk);
end