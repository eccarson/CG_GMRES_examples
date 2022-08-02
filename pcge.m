
% Run the CG method to solve Ax=b in "exact" arithmetic
% Note: Requires Advanpix library

%Input:
%A: square, sparse matrix with dimension n
%b: right hand side of system to solve, Ax=b; vector of dimension n
%P: preconditioner to be applied (left preconditioning)
%tol: convergence criteria for computed residual 2-norm
%options: should be a struct that contains two values:
%   options.xlim: x axis limit for plotting
%   options.truesol: true solution (used for computing the error)

%Output:
%results struct stores:
%r_exact_norm: 2-norm of true residual computed in each iteration
%(results.r_exact_norm)
%r_comp_norm: 2-norm of computed residual computed in each iteration
%(results.r_comp_norm)
%error_A_norm: A-norm of the error
%(results.error_A_norm)
%x: approximate solution computed in each iteration
%(results.x)
%r: residual computed in each iteration
%(results.r)
%ritzvals: ritz values of the Jacobi matrix T
%(results.ritzvals)

function results = pcge(A, b, P, tol, options)

% 200 digits of precision to simulate exact arithmetic
mp.Digits(200);

A = mp(A); P = mp(P); b = mp(b);

% Size of matrix
N = size(A,1);

x0 = mp(zeros(N,1));

% Set initial values for vectors
r0 = b - A*x0;
p0 = r0;
x(:,1)  = x0;
r(:,1)  = r0;
p(:,1)  = p0;
z(:,1) = P\r0;

% Set total number of iterations to 0
its = 0;

% Initialize initial true and computed residuals
results.r_exact_norm(1) = norm(b-A*x0);
results.r_comp_norm(1) = norm(r0);
results.error_A_norm(1) = 1;
results.x=x0;


% Begin the iterations
while its < options.xlim
    
    %     %Break out of the loop if we have converged
    %     if(results.r_comp_norm(its+1) <= tol)
    %         break;
    %     end
    
    % increase iteration count
    its = its + 1;
    
    % Compute scalar alpha
    alpha(its) = r(:,its)'*z(:,its)/(p(:,its)'*A*p(:,its));
    
    % Update x coordinate vector
    x(:,its+1) = x(:,its) + alpha(its)*p(:,its);
    
    % Update r coordinate vector
    r(:,its+1) = r(:,its) - alpha(its)*A*p(:,its);
    
    z(:,its+1) = P\r(:,its+1);
    
    % Compute scalar beta
    beta(its) = (r(:,its+1)'*z(:,its+1))/ (r(:,its)'*z(:,its));
    
    % Update p coordinate vector
    p(:,its+1) = z(:,its+1) + beta(its)*p(:,its);
    
    % Compute and store true residual norm
    results.r_exact_norm(its+1) = norm(b-A*x(:,its+1));
    
    % Compute and store computed residual norm
    results.r_comp_norm(its+1) = norm(r(:,its+1));
    
    err = mp(x(:,its+1) - options.truesol,34);
    num = sqrt(double(mp(err',34)*mp(A,34)*mp(err,34)));
    denom = sqrt(double(mp(options.truesol',34)*mp(A,34)*mp(options.truesol,34)));
    results.error_A_norm(its+1) = double(num/denom);
    results.r = r;
    
    % Find Ritz values
    L = diag(1./sqrt(alpha(1:its)),0) + diag(sqrt(beta(1:its-1)./alpha(1:its-1)),1);
    T = L'*L;   
    results.ritzvals = eig(full(T));
    
    
    % Store current solution
    results.x = x(:,its+1);
    
end






% % Run the CG method to solve Ax=b in "exact" arithmetic
% % Note: Requires Advanpix library
% 
% %Input:
% %A: square, sparse matrix with dimension n
% %b: right hand side of system to solve, Ax=b; vector of dimension n
% %P: preconditioner to be applied (left preconditioning)
% %tol: convergence criteria for computed residual 2-norm
% %options: should be a struct that contains two values:
% %   options.xlim: x axis limit for plotting
% %   options.truesol: true solution (used for computing the error)
% 
% %Output:
% %results struct stores:
% %r_exact_norm: 2-norm of true residual computed in each iteration
% %(results.r_exact_norm)
% %r_comp_norm: 2-norm of computed residual computed in each iteration
% %(results.r_comp_norm)
% %error_A_norm: A-norm of the error
% %(results.error_A_norm)
% %x: approximate solution computed in each iteration
% %(results.x)
% %r: residual computed in each iteration
% %(results.r)
% %ritzvals: ritz values of the Jacobi matrix T
% %(results.ritzvals)
% 
% function results = pcge(A, b, P, tol, options)
% 
% 
% % Size of matrix
% N = size(A,1);
% x0 = zeros(N,1);
% 
% % Apply preconditioner
% A = P\A;
% b = P\b;
% 
% 
% % Run Lanczos with double reorthogonalization
% [V,T,bet] = tridiag_lan_double(A,b-A*x0,options.xlim,2);
% results.ritzvals = eig(full(T));
% 
% % Compute solutions and residuals
% [Xk,Rk,rescomp] = comp_solutions_lanczos(T,V,bet);
% results.r = Rk;
% results.x = Xk(:,end);
% results.r_comp_norm = rescomp;
% 
% % Compute true residuals
% [Rkt,restrue] = comp_trueres(A,b,Xk);
% results.r_exact_norm = restrue;
% 
% % Compute errors
% [err] = comp_error(A,options.truesol,Xk);
% err = err./sqrt(options.truesol'*A*options.truesol);
% results.error_A_norm = err;
% 


