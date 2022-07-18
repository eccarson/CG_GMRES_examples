% This function generates convergence plots for the example:
% "Worst-case CG"

% Set 100 digits to simulate "exact" arithmetic
mp.Digits(100);

% Generate diagonal matrix and right-hand side
n = 48;
l1 = 1;
ln = 5;
A = mp(strakosmatrix(n, l1, ln, 1),100);
b = ones(n,1)./sqrt(n);

% Compute kappa-based upper bound on convergence rate
kappa = ln/l1;
factor = (sqrt(kappa)-1)/(sqrt(kappa)+1);
for k = 0:n
    kappa_bd(k+1) = 2*factor^k;
end

% Compute true solution and set options for cg run
options.truesol = double(mp(A)\mp(b));
tol = 1e-16; options.xlim = n;

% Run exact CG on this problem
results = pcge(A, b, eye(n), tol, options);
Anorm = results.error_A_norm;

%Plot A-norm of the error
figure(1), clf
figure(1), semilogy([0:length(Anorm)-1],Anorm,'r--', 'LineWidth',2), hold on
figure(1), semilogy([0:length(kappa_bd)-1],kappa_bd,'k-', 'LineWidth',2)

axis([0 n 1e-12 10])
ax = gca;
set(ax,'FontSize',16)
print -depsc worst1.eps


% Generate second problem; same as first, but with larger ln. 
n = 48;
l1 = 1;
ln = 100;
A = mp(strakosmatrix(n, l1, ln, 1),100);
b = ones(n,1)./sqrt(n);

% Compute kappa-based upper bound on convergence rate
kappa = ln/l1;
factor = (sqrt(kappa)-1)/(sqrt(kappa)+1);
for k = 0:n
    kappa_bd(k+1) = 2*factor^k;
end

% Compute true solution and set options for cg run
options.truesol = double(mp(A)\mp(b));
tol = 1e-16; options.xlim = n;

% Run exact CG on this problem
results = pcge(A, b, eye(n), tol, options);
Anorm = results.error_A_norm;

%Plot A-norm of the error
figure(2), clf
figure(2), semilogy([0:length(Anorm)-1],Anorm,'r--', 'LineWidth',2), hold on
figure(2), semilogy([0:length(kappa_bd)-1],kappa_bd,'k-', 'LineWidth',2)
axis([0 n 1e-12 10])
ax = gca;
set(ax,'FontSize',16)
print -depsc worst2.eps


% Generate third problem; same as first, but with different eigenvalue
% distribution
n = 48;
l1 = 1;
ln = 5;
A = mp(strakosmatrix(n, l1, ln, .1),100);
b = ones(n,1)./sqrt(n);

% Compute kappa-based upper bound on convergence rate
kappa = ln/l1;
factor = (sqrt(kappa)-1)/(sqrt(kappa)+1);
for k = 0:n
    kappa_bd(k+1) = 2*factor^k;
end

% Compute true solution and set options for cg run
options.truesol = double(mp(A)\mp(b));
tol = 1e-16; options.xlim = n;

% Run exact CG on this problem
results = pcge(A, b, eye(n), tol, options);
Anorm = results.error_A_norm;

%Plot A-norm of the error
figure(3), clf
figure(3), semilogy([0:length(Anorm)-1],Anorm,'r--', 'LineWidth',2), hold on
figure(3), semilogy([0:length(kappa_bd)-1],kappa_bd,'k-', 'LineWidth',2)

axis([0 n 1e-12 10])
ax = gca;
set(ax,'FontSize',16)
print -depsc worst3.eps


% Generate fourth problem; same as first, but with different RHS
n = 48;
l1 = 1;
ln = 5;
A = mp(strakosmatrix(n, l1, ln, 1),100);
b = [1;1;zeros(n-2,1)];
b = b./norm(b);

% Compute kappa-based upper bound on convergence rate
kappa = ln/l1;
factor = (sqrt(kappa)-1)/(sqrt(kappa)+1);
for k = 0:n
    kappa_bd(k+1) = 2*factor^k;
end

% Compute true solution and set options for cg run
options.truesol = double(mp(A)\mp(b));
tol = 1e-16; options.xlim = 2;

% Run exact CG on this problem
results = pcge(A, b, eye(n), tol, options);
Anorm = results.error_A_norm;

%Plot A-norm of the error
figure(4), clf
figure(4), semilogy([0:length(Anorm)-1],Anorm,'r--', 'LineWidth',2), hold on
figure(4), semilogy([0:length(kappa_bd)-1],kappa_bd,'k-', 'LineWidth',2)
axis([0 n 1e-12 10])
ax = gca;
set(ax,'FontSize',16)
print -depsc worst4.eps
