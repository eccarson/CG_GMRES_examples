% This function generates the spectral density plots for the example:
% "Mathematical behavior of CG for problems with clustered eigenvalues"

function spectraldensityclustexact

%Generate problems with 10 eigenvalues
n = 10;
l1 = .1;
ln = 1e3;
rho = 0.6;

% Construct 3 different eigenvalue distributions
lam1 = diag(strakosmatrix(n, l1, ln, rho));
lam2 = diag(strakosmatrix2(n, l1, ln, rho));
lam3 = diag(strakosmatrix(n, l1, ln, 1));


% Create clusters
nclust = 10;
intv = 1e-12;
lam1c = zeros(n*nclust,1); lam2c = zeros(n*nclust,1); lam3c = zeros(n*nclust,1);
for i = 1:n
    for j = 1:nclust/2
        lam1c((i-1)*nclust+2*j-1) = lam1(i)+j*intv;
        lam2c((i-1)*nclust+2*j-1) = lam2(i)+j*intv;
        lam3c((i-1)*nclust+2*j-1) = lam3(i)+j*intv;
        lam1c((i-1)*nclust+2*j) = lam1(i)-j*intv;
        lam2c((i-1)*nclust+2*j) = lam2(i)-j*intv;
        lam3c((i-1)*nclust+2*j) = lam3(i)-j*intv;
    end
end

lam1c = sort(lam1c);
lam2c = sort(lam2c);
lam3c = sort(lam3c);

% Construct diagonal matrices and right-hand sides
lam1c = diag(lam1c); lam2c = diag(lam2c); lam3c = diag(lam3c);
b = ones(size(lam1c,1),1);
b = b./norm(b);

% CG tolerance
tol = 1e-14;

% 200 digits of precision to simulate exact arithmetic; only used here for
% computing true solutions
mp.Digits(200);


% Compute true solutions
optionslam1.truesol = double(mp(lam1c)\mp(b));
optionslam2.truesol = double(mp(lam2c)\mp(b));
optionslam3.truesol = double(mp(lam3c)\mp(b));

% axis limit parameters for plotting
optionslam1.xlim=5; 
optionslam2.xlim=5;
optionslam3.xlim=5;

% Run exact CG for the systems; 5 iterations
resultslam1d5 = pcge(double(lam1c), double(b), eye(size(lam1c,1)), tol, optionslam1);
resultslam2d5 = pcge(double(lam2c), double(b), eye(size(lam2c,1)), tol, optionslam2);
resultslam3d5 = pcge(double(lam3c), double(b), eye(size(lam3c,1)), tol, optionslam3);

optionslam1.xlim=10; 
optionslam2.xlim=10; 
optionslam3.xlim=10;

% Run exact CG for the systems; 10 iterations
resultslam1d10 = pcge(double(lam1c), double(b), eye(size(lam1c,1)), tol, optionslam1);
resultslam2d10 = pcge(double(lam2c), double(b), eye(size(lam2c,1)), tol, optionslam2);
resultslam3d10 = pcge(double(lam3c), double(b), eye(size(lam3c,1)), tol, optionslam3);

optionslam1.xlim=15;

% Run exact CG for the systems; 15 iterations
resultslam1d15 = pcge(double(lam1c), double(b), eye(size(lam1c,1)), tol, optionslam1);


% Plot spectral densities
figure(1)
subplot(3,1,2)


stairs(diag(lam1c),(1/(n*nclust)):(1/(n*nclust)):1,  'r-', 'LineWidth',1);
hold on;
stairs(resultslam1d5.ritzvals, (1/5):(1/5):1, ":", 'LineWidth',1, 'Color', [0,0,0]);
stairs(resultslam1d10.ritzvals,(1/10):(1/10):1,  "-.", 'LineWidth',1,'Color', [0,0,0]);
stairs(resultslam1d15.ritzvals,(1/15):(1/15):1,  "--", 'LineWidth',1, 'Color', [0,0,0]);
legend('eigenvalues', 'Ritz values k = 5', 'Ritz values k = 10', 'Ritz values k = 15','Location','Southeast','Interpreter','latex');
title('acc. to the left','Interpreter','latex');
axis([0, 1000, 0, 1])
set(gca, 'fontsize', 14);



subplot(3,1,1)
stairs( diag(lam2c),(1/(n*nclust)):(1/(n*nclust)):1, 'b-', 'LineWidth',1);
hold on;
stairs(resultslam2d5.ritzvals, (1/5):(1/5):1, ":", 'LineWidth',1, 'Color', [0,0,0]);
stairs(resultslam2d10.ritzvals,(1/10):(1/10):1,  "-.", 'LineWidth',1, 'Color', [0,0,0]);
legend('eigenvalues', 'Ritz values k = 5','Ritz values k = 10','Location','Northwest','Interpreter','latex');
title('acc. to the right','Interpreter','latex');
axis([0, 1000, 0, 1])
set(gca, 'fontsize', 14);



subplot(3,1,3)
stairs(diag(lam3c),(1/(n*nclust)):(1/(n*nclust)):1,  'g-', 'LineWidth',1);
hold on;
stairs(resultslam3d5.ritzvals,(1/5):(1/5):1,  ":", 'LineWidth',1, 'Color', [0,0,0]);
stairs(resultslam3d10.ritzvals,(1/10):(1/10):1,  "-.", 'LineWidth',1, 'Color', [0,0,0]);
legend('eigenvalues', 'Ritz values k = 5','Ritz values k = 10','Ritz values k = 30', 'Location','Northwest','Interpreter','latex');
title('equally spaced','Interpreter','latex');
axis([0, 1000, 0, 1])
set(gca, 'fontsize', 14);
% print -depsc specdens.eps



