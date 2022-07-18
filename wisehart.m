% This script generates Wisehart eigenvalues and convergence plot for the example:
% "Worst-case CG"

% Joerg Liesen, 12.02.2020

format compact

n = 500; m = 100;
lam = (m-1)/n;
cucker_val = exp(-m*(1-lam)/91);
i = sqrt(-1);

tol = 1e-10; MAXIT = m/2; x0 = zeros(m,1);

f1 =figure;
figure(f1), clf, hold on
f2 = figure;
figure(f2), clf

% For 100 trials, generate random Wisehart matrices and plot their spectra
% and CG convergence trajectory
for k = 1:100
    A = randn(n,m); A = A'*A;
    condA(k) = cond(A);
    f = eig(A);
    figure(f1), plot(real(f),imag(f+k*i),'.')
    x = randn(m,1); b = A*x;
    [Anorm] = my_cg(A,b,x0,x,tol,MAXIT);
    figure(f2), semilogy([0:length(Anorm)-1],Anorm,'--'), hold on
end

figure(1)
ax = gca;
set(ax,'FontSize',16)

mean_condA = mean(condA);

% Compute and plot kappa-based upper bound
kappa = (sqrt(mean_condA)-1)/(sqrt(mean_condA)+1);
kappa_bd(1) = 1;
for k=1:(m/2)-1
    kappa_bd(k+1) = 2*(1/(kappa^k+kappa^(-k)));
end
figure(f2), semilogy([0:length(kappa_bd)-1],kappa_bd,'k','LineWidth',2)
axis([0 length(kappa_bd)-1 1e-16 1])
ax = gca;
set(ax,'FontSize',16)

% Save figures
figure(f1), print -depsc evals_SPD_rand.eps
figure(f2), print -depsc CG_SPD_rand.eps

