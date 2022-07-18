% This script plots convergence trajectories for the example:
% "CG with "average" initial residuals"

clear all

% Set digits to simulate exact arithmetic
mp.Digits(100);

% Construct two diagonal matrices with different eigenvalue distributions
n = 48;
l1 = .1;
ln = 1e3;
A = mp(strakosmatrix(n, l1, ln, .6),100);
A2 = mp(strakosmatrix(n, l1, ln, 1),100);

figure()

% Set options for CG run
tol = 1e-16; options.xlim = n; options2.xlim = n; 
x0 = zeros(n,1);

% Seed random number generator for reproducibility
rng(1);

% Perform 100 trials
for k = 1:100
    
    % Generate random RHS
    b = randn(n,1);
    
    % Compute true solution for first system
    options.truesol = (mp(A)\mp(b));
    
    % Run exact CG
    results = pcge_x0(A, b, eye(n), x0, tol, options);
    Anorm1 = results.error_A_norm;
    
    % Compute true solution for second system
    options2.truesol = (mp(A2)\mp(b));
    
    % Run exact CG
    results2 = pcge_x0(A2, b, eye(n), x0, tol, options2);
    Anorm2 = results2.error_A_norm;
    
    % Plot error trajectories
    semilogy([0:length(Anorm1)-1],Anorm1,'r:'), hold on
    semilogy([0:length(Anorm2)-1],Anorm2,'b:'), 
end

axis([0 n 1e-12 10])
ax = gca;
set(ax,'FontSize',16)
print -depsc randx0.eps
