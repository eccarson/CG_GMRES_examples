% This function generates the spectral density plots for the example:
% "Computational behavior of CG for different eigenvalue distributions"

function spectraldensity

n = 30;
l1 = .1;
ln = 1e3;
rho = .6;

% Construct 3 different eigenvalue distributions
lam1 = diag(strakosmatrix(n, l1, ln, rho));
lam2 = diag(strakosmatrix2(n, l1, ln, rho));
lam3 = diag(strakosmatrix(n, l1, ln, 1));

% 200 digits of precision to simulate exact arithmetic; only used for
% computing true solutions
mp.Digits(200);

% Construct diagonal matrices and right-hand sides
lam1 = diag(lam1); lam2 = diag(lam2); lam3 = diag(lam3);
b = ones(size(lam1,1),1);
b = b./norm(b);

tol = 1e-14;

% Compute true solutions
optionslam1.truesol = double(mp(lam1)\mp(b));
optionslam2.truesol = double(mp(lam2)\mp(b));
optionslam3.truesol = double(mp(lam3)\mp(b));

% axis limit parameters for plotting
optionslam1.xlim=10; 
optionslam2.xlim=10;
optionslam3.xlim=10;

% Run double precision CG for the systems; 10 iterations
resultslam1d10 = pcgd(double(lam1), double(b), eye(size(lam1,1)), tol, optionslam1);
resultslam2d10 = pcgd(double(lam2), double(b), eye(size(lam2,1)), tol, optionslam2);
resultslam3d10 = pcgd(double(lam3), double(b), eye(size(lam3,1)), tol, optionslam3);

optionslam1.xlim=20; 
optionslam3.xlim=20;

% Run double precision CG for the systems; 20 iterations
resultslam1d20 = pcgd(double(lam1), double(b), eye(size(lam1,1)), tol, optionslam1);
resultslam3d20 = pcgd(double(lam3), double(b), eye(size(lam3,1)), tol, optionslam3);

optionslam1.xlim=30; 
optionslam3.xlim=30;

% Run double precision CG for the systems; 30 iterations
resultslam1d30 = pcgd(double(lam1), double(b), eye(size(lam1,1)), tol, optionslam1);
resultslam3d30 = pcgd(double(lam3), double(b), eye(size(lam3,1)), tol, optionslam3);


% Plot spectral densities
figure(1)
subplot(3,1,2)

stairs(diag(lam1),0:(1/(n-1)):1,  'r-', 'LineWidth',1);
hold on;
stairs(resultslam1d10.ritzvals, 0:(1/9):1, "--", 'LineWidth',1, 'Color', [.8, .8, .8]);
stairs(resultslam1d20.ritzvals,0:(1/19):1,  "--", 'LineWidth',1,'Color', [.4, .4, .4]);
stairs(resultslam1d30.ritzvals,0:(1/29):1,  "--", 'LineWidth',1, 'Color', [0, 0, 0]);
%stairs(resultslam1d50.ritzvals,0:(1/49):1,  "--", 'LineWidth',1, 'Color', [0, 0,0]);
legend('eigenvalues', 'Ritz values k = 10', 'Ritz values k = 20', 'Ritz values k = 30','Location','Southeast','Interpreter','latex');
title('acc. to the left','Interpreter','latex');
set(gca, 'fontsize', 14);

subplot(3,1,1)
stairs( diag(lam2),0:(1/(n-1)):1, 'b-', 'LineWidth',1);
hold on;
stairs(resultslam2d10.ritzvals,0:(1/9):1,  "--", 'LineWidth',1, 'Color', [.8, .8, .8]);
legend('eigenvalues', 'Ritz values k = 10','Location','Northwest','Interpreter','latex');
title('acc. to the right','Interpreter','latex');
set(gca, 'fontsize', 14);

subplot(3,1,3)
stairs(diag(lam3),0:(1/(n-1)):1,  'g-', 'LineWidth',1);
hold on;
stairs(resultslam3d10.ritzvals,0:(1/9):1,  "--", 'LineWidth',1, 'Color', [.8, .8, .8]);
stairs(resultslam3d20.ritzvals,0:(1/19):1,  "--", 'LineWidth',1, 'Color', [.4, .4, .4]);
stairs(resultslam3d30.ritzvals,0:(1/29):1,  "--", 'LineWidth',1, 'Color', [0,0,0]);
legend('eigenvalues', 'Ritz values k = 10','Ritz values k = 20','Ritz values k = 30', 'Location','Northwest','Interpreter','latex');
title('equally spaced','Interpreter','latex');
set(gca, 'fontsize', 14);
% print -depsc specdens.eps





