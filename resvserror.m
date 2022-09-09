% This script creates convergence plots for the example:
% "Residual Versus Error and Stopping Criteria"
%
% Note: Construction from:
% [Meurant, Gérard. "On prescribing the convergence behavior of the 
% conjugate gradient algorithm." Numerical Algorithms 84.4 (2020): 
% 1353-1380].

% Construct example where residual doesn't decrease but error does
n = 20;
resnorms = [1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2];
errors(1) = 1;
for i = 2:20
    errors(i) = errors(i-1)*.4;
end

for i = 1:n
    nu(i) = 1/resnorms(i);
    sigma(i) = errors(i)^2/(resnorms(i)*resnorms(1)); %errors(i)^2*nu(i)/(resnorms(1)^2);
end

for i = 1:n
   for j = 1:i
       if j == 1
           L(i,j) = sigma(i);
       else
           L(i,j) = sigma(i)*nu(j);
       end
   end
end

Tinv = L + tril(L,-1)';
T = inv(mp(Tinv,200));
A = T; 
b = [1; zeros(n-1,1)];


% Run exact CG for this system
options.xlim = n;
options.truesol = double(mp(A,200)\mp(b,200));
results = pcge(A, b, speye(n), 1e-17, options);

% Plot convergence of residual and error
figure()
semilogy(0:numel(results.r_exact_norm)-1, results.r_exact_norm, 'b-', 'LineWidth',1);
hold on;
semilogy(0:numel(results.error_A_norm)-1, results.error_A_norm, 'r--', 'LineWidth',1);
axis([0,19,1e-8,1e1])
% legend('residual, $\Vert b-Ax_k\Vert_2$', 'error, $\Vert x-x_k\Vert_A$', 'Location','Southwest','Interpreter','latex');
set(gca, 'fontsize', 16);
print -depsc resvserror1.eps


% Construct example where error doesn't decrease fast but residual does
n = 20;
resnorms(1) = 1; 
errors(1) = 1;
for i = 2:20
    errors(i) = errors(i-1)*.999;
    resnorms(i) = resnorms(i-1)*.4;
end

for i = 1:n
    nu(i) = 1/resnorms(i);
    sigma(i) = errors(i)^2*nu(i)/(resnorms(1)^2);
end

for i = 1:n
   for j = 1:i
       if j == 1
           L(i,j) = sigma(i);
       else
           L(i,j) = sigma(i)*nu(j);
       end
   end
end

Tinv = L + tril(L,-1)';
T = inv(mp(Tinv,200));
A = T; 
b = [1; zeros(n-1,1)];

% Run exact CG for this system
options.xlim = n;
options.truesol = double(mp(A,200)\mp(b,200));
results = pcge(A, b, speye(n), 1e-17, options);

% Plot convergence of residual and error
figure()
semilogy(0:numel(results.r_exact_norm)-1, results.r_exact_norm, 'b-', 'LineWidth',1);
hold on;
semilogy(0:numel(results.error_A_norm)-1, results.error_A_norm, 'r--', 'LineWidth',1);
axis([0,19,1e-8,1e1])
% legend('residual, $\Vert b-Ax_k\Vert_2$', 'error, $\Vert x-x_k\Vert_A$', 'Location','Southwest','Interpreter','latex');
set(gca, 'fontsize', 16);
print -depsc resvserror2.eps