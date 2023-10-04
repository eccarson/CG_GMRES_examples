%Compute true solutions at each iteration given matrix A, RHS b,
%and matrix X whose columns are the approximate solution vectors at 
%each iteration
%Returns matrix Rk whose columns are the true solution at each iteration 
%and vector res which stores the 2-norms of these columns; computations all
%in simulated exact arithmetic
%
% Note: Requires Advanpix library
%
function [Rk,res] = comp_trueres_exact(A,b,X)

% 200 digits of precision to simulate exact arithmetic
mp.Digits(200);


A = mp(A);
b = mp(b);
X = mp(X);

for i = 1:size(X,2)-1
   Rk(:,i) = b-A*X(:,i);
   res(i) = norm(Rk(:,i));
end
