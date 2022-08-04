% This function, run with rho=0.6, generates the convergence plots for the examples:
% 1) "Mathematical behavior of CG for problems with clustered eigenvalues"

function clusterexp(rho)

%Generate problems with 10 eigenvalues
n = 10;
l1 = .1;
ln = 1e3;

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



% Construct diagonal matrices and right-hand sides
lam1c = diag(lam1c); lam2c = diag(lam2c); lam3c = diag(lam3c);
b = ones(size(lam1c,1),1);
b = b./norm(b);

% CG tolerance
tol = 1e-14;

% axis limit parameters for plotting
optionslam1.xlim=80; 
optionslam2.xlim=80;
optionslam3.xlim=60;

mp.Digits(100);

% Run exact CG for the 3 systems
optionslam1.truesol = double(mp(lam1c)\mp(b));
resultslam1 = pcge(lam1c, b, eye(size(lam1c,1)), tol, optionslam1);

optionslam2.truesol = double(mp(lam2c)\mp(b));
resultslam2 = pcge(lam2c, b, eye(size(lam2c,1)), tol, optionslam2);

optionslam3.truesol = double(mp(lam3c)\mp(b));
resultslam3 = pcge(lam3c, b, eye(size(lam3c,1)), tol, optionslam3);

% cut off extra points for plotting
i = find(resultslam1.error_A_norm < tol,1,'first');
if(~isempty(i))
    resultslam1.error_A_norm = resultslam1.error_A_norm(1:i);
end

i = find(resultslam2.error_A_norm < tol,1,'first');
if(~isempty(i))
    resultslam2.error_A_norm = resultslam2.error_A_norm(1:i);
end

i = find(resultslam3.error_A_norm < 7e-16,1,'first');
if(~isempty(i))
    resultslam3.error_A_norm = resultslam3.error_A_norm(1:i);
end


% Plot convergence for exact CG
figure()
semilogy(0:numel(resultslam2.error_A_norm)-1,resultslam2.error_A_norm,'b--','LineWidth',2);
hold on;
semilogy(0:numel(resultslam1.error_A_norm)-1,resultslam1.error_A_norm,'r--','LineWidth',2);
semilogy(0:numel(resultslam3.error_A_norm)-1,resultslam3.error_A_norm,'g--','LineWidth',2);

axis([0,20,4e-11,10])
legend('acc. to the right','acc. to the left',  'equally spaced','Interpreter','latex');
set(gca,'FontSize',16)
print -depsc exactcgclust.eps





