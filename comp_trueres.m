%Compute true solutions at each iteration given matrix A, RHS b,
%and matrix X whose columns are the approximate solution vectors at 
%each iteration
%Returns matrix Rk whose columns are the true solution at each iteration 
%and vector res which stores the 2-norms of these columns; computations all
%in double precision
function [Rk,res] = comp_trueres(A,b,X)

A = double(A);
b = double(b);
X = double(X);

for i = 1:size(X,2)-1
   Rk(:,i) = b-A*X(:,i);
   res(i) = norm(Rk(:,i));
end
