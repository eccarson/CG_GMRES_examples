%Run Lanczos tridiagonalization in double precision
%with matrix A, vector b, for k iterations
%reo = number of reorthogonalization steps
%Outputs V and T matrices in Lanczos recurrence
function [V,T,bet] = tridiag_lan_double(A,b,k,reo)

n = size(A,2);
T       = zeros(k,k);
V       = zeros(n,k);

beta    = 0.0; 
w      = b;
q1     = b;

bet = norm(b);


for j = 1 : k
    
    q = w/sqrt(w'*w); 
    V(:,j) = q;
    w = A*q - sqrt(beta)*q1;
    alpha = q'*w;
    w = w - alpha*q;
    beta = w'*w;
    q1 = q;
%
    for count = 1 : reo
        for k = 1 : j
            w = w - (V(:,k)'*w)*V(:,k);
        end
    end
%
    T(j,j)    = alpha;
    T(j,j+1)  = sqrt(beta);
    T(j+1,j)  = sqrt(beta);  
end
q = w/sqrt(w'*w); 
V(:,k+1) = q;

