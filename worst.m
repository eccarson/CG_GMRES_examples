% This function generates convergence plots for the example:
% "Worst-case CG"

% Set 100 digits to simulate "exact" arithmetic
mp.Digits(200);

% Generate diagonal matrix and right-hand side
n = 48;
l1 = 1;
ln = 5;
A = mp(strakosmatrix(n, l1, ln, 1),100);
b = ones(n,1)./sqrt(n);

% Compute kappa-based upper bound on convergence rate and minmax bound
kappa = ln/l1;
lam = double(eig(full(A)));
out(1)=1; kappa_bd(1)=2;
factor = (sqrt(kappa)-1)/(sqrt(kappa)+1);
for k = 1:n-1
    [blk,Avec,C,bb,X0,y0,Z0,objval,p] = cheby0(lam,k,1);
    out(k+1)=objval;
    kappa_bd(k+1) = 2*factor^k;
end

% Compute true solution and set options for cg run
options.truesol = double(mp(A)\mp(b));
tol = 1e-16; options.xlim = n;

% Run exact CG on this problem
results = pcge(A, b, eye(n), tol, options);
Anorm = results.error_A_norm;

%Plot A-norm of the error
f1 = figure;
figure(f1), semilogy([0:length(Anorm)-1],Anorm,'b-', 'LineWidth',2), hold on
figure(f1), semilogy([0:length(out)-1],out,'k:','LineWidth',2)
figure(f1), semilogy([0:length(kappa_bd)-1],kappa_bd,'k--', 'LineWidth',2)

title('case 1','Interpreter','latex');
axis([0 40 1e-8 10])
ax = gca;
set(ax,'FontSize',16)
print -depsc worst1.eps


% Generate second problem; same as first, but with larger ln. 
n = 48;
l1 = 1;
ln = 100;
A = mp(strakosmatrix(n, l1, ln, 1),100);
b = ones(n,1)./sqrt(n);

% Compute kappa-based upper bound on convergence rate and minmax bound
kappa = ln/l1;
lam = double(eig(full(A)));
out(1)=1; kappa_bd(1)=2;
factor = (sqrt(kappa)-1)/(sqrt(kappa)+1);
for k = 1:n-1
    [blk,Avec,C,bb,X0,y0,Z0,objval,p] = cheby0(lam,k,1);
    out(k+1)=objval;
    kappa_bd(k+1) = 2*factor^k;
end

% Compute true solution and set options for cg run
options.truesol = double(mp(A)\mp(b));
tol = 1e-16; options.xlim = n;

% Run exact CG on this problem
results = pcge(A, b, eye(n), tol, options);
Anorm = results.error_A_norm;

%Plot A-norm of the error
f2 = figure;
figure(f2), semilogy([0:length(Anorm)-1],Anorm,'b-', 'LineWidth',2), hold on
figure(f2), semilogy([0:length(out)-1],out,'k:','LineWidth',2)
figure(f2), semilogy([0:length(kappa_bd)-1],kappa_bd,'k--', 'LineWidth',2)
title('case 2','Interpreter','latex');
axis([0 40 1e-8 10])
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

% Compute kappa-based upper bound on convergence rate and minmax bound
kappa = ln/l1;
lam = double(eig(full(A)));
out(1)=1; kappa_bd(1)=2;
factor = (sqrt(kappa)-1)/(sqrt(kappa)+1);
for k = 1:n-1
    [blk,Avec,C,bb,X0,y0,Z0,objval,p] = cheby0(lam,k,1);
    out(k+1)=objval;
    kappa_bd(k+1) = 2*factor^k;
end

% Compute true solution and set options for cg run
options.truesol = double(mp(A)\mp(b));
tol = 1e-16; options.xlim = n;

% Run exact CG on this problem
results = pcge(A, b, eye(n), tol, options);
Anorm = results.error_A_norm;

%Plot A-norm of the error
f3 = figure;
figure(f3), semilogy([0:length(Anorm)-1],Anorm,'b-', 'LineWidth',2), hold on
figure(f3), semilogy([0:length(out)-1],out,'k:','LineWidth',2)
figure(f3), semilogy([0:length(kappa_bd)-1],kappa_bd,'k--', 'LineWidth',2)
title('case 3','Interpreter','latex');
axis([0 40 1e-8 10])
ax = gca;
set(ax,'FontSize',16)
print -depsc worst3.eps


% Generate fourth problem; same as first, but with different RHS
n = 48;
l1 = 1;
ln = 5;
A = mp(strakosmatrix(n, l1, ln, 1),100);
b = [1;1e-12.*ones(n-2,1);1];
b = b./norm(b);

% Compute kappa-based upper bound on convergence rate and minmax bound
kappa = ln/l1;
lam = double(eig(full(A)));
out(1)=1; kappa_bd(1)=2;
factor = (sqrt(kappa)-1)/(sqrt(kappa)+1);
for k = 1:n-1
    [blk,Avec,C,bb,X0,y0,Z0,objval,p] = cheby0(lam,k,1);
    out(k+1)=objval;
    kappa_bd(k+1) = 2*factor^k;
end
% Compute true solution and set options for cg run
options.truesol = double(mp(A)\mp(b));
tol = 1e-16; options.xlim = 10;

% Run exact CG on this problem
results = pcge(A, b, eye(n), tol, options);
Anorm = results.error_A_norm;


%Plot A-norm of the error
f4 = figure;
figure(f4), semilogy([0:length(Anorm)-1],Anorm,'b-', 'LineWidth',2), hold on
figure(f4), semilogy([0:length(out)-1],out,'k:','LineWidth',2)
figure(f4), semilogy([0:length(kappa_bd)-1],kappa_bd,'k--', 'LineWidth',2)
title('case 4','Interpreter','latex');
axis([0 40 1e-8 10])
ax = gca;
set(ax,'FontSize',16)
print -depsc worst4.eps
