%Compute (in simulated exact arithmetic) the CGL solutions 
%and "computed" residuals 
%in each iteration given the V and T matrices 
%output from Lanczos and parameter beta
%
% Note: Requires Advanpix library
%
function [Xk,Rk,res] = comp_solutions_lanczos_exact(T,V,beta)


% 200 digits of precision to simulate exact arithmetic
mp.Digits(200);

T = mp(T);
V = mp(V);
beta = mp(beta);

rhs = mp(zeros(length(T),1)); rhs(1) = beta;

Xk = mp(zeros(size(V,1),size(T,2)));%to store the solutions
Rk = mp(zeros(size(V,1),size(T,2)));%to store the residuals

for k = 1:size(T,2)-1 %compute solutions exact

    yk = T(1:k,1:k)\rhs(1:k);
    xk = V(:,1:k)*yk;
    e1 = mp([1;zeros(k,1)]);
    rk = V(:,1:k+1)*(beta*e1-T(1:k+1,1:k)*yk);
        
    Xk(:,k) = xk;%store solutions
    Rk(:,k) = rk;%store residuals
    res(k) = norm(rk);
    
end