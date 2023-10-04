% Run Lanczos tridiagonalization in simulated exact arithmetic
% with matrix A, vector b, for k iterations
% Outputs V and T matrices in Lanczos recurrence
%
% Note: Requires Advanpix library
%
function [V,T,bet] = tridiag_lan_exact(A,b,k)

% 200 digits of precision to simulate exact arithmetic
mp.Digits(200);

A = mp(A); b = mp(b);

n = size(A,2);
T       = mp(zeros(k,k));
V       = mp(zeros(n,k));

beta    = mp(0.0); 
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
    for k = 1 : j
            w = w - (V(:,k)'*w)*V(:,k);
    end
%
    T(j,j)    = alpha;
    T(j,j+1)  = sqrt(beta);
    T(j+1,j)  = sqrt(beta);  
end
q = w/sqrt(w'*w); 
V(:,k+1) = q;

