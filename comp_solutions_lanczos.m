%Compute the CGL solutions and "computed" residuals 
%in each iteration given the V and T matrices 
%output from Lanczos and parameter beta
function [Xk,Rk,res] = comp_solutions_lanczos(T,V,beta)

T = double(T);
V = double(V);
beta = double(beta);

rhs = zeros(length(T),1); rhs(1) = beta;

Xk = zeros(size(V,1),size(T,2));%to store the solutions
Rk = zeros(size(V,1),size(T,2));%to store the residuals

for k = 1:size(T,2)-1 %compute solutions exact

    yk = T(1:k,1:k)\rhs(1:k);
    xk = V(:,1:k)*yk;
    e1 = [1;zeros(k,1)];
    rk = V(:,1:k+1)*(beta*e1-T(1:k+1,1:k)*yk);
        
    Xk(:,k) = xk;%store solutions
    Rk(:,k) = rk;%store residuals
    res(k) = norm(rk);
    
end