% This script generates convergence plot for the example:
% "Computational behavior of different CG implementations"

% Set digits to simulate exact arithmetic
mp.Digits(100);

% Construct digonal matrix and right-hand side
n = 48;
l1 = .1;
ln = 1e3;
lam1 = strakosmatrix(n, l1, ln, .25);
b1 = ones(n,1);
b1 = b1./norm(b1);

% Set parameters for CG run
tol = 6e-17;
optionslam1.xlim=80; 
nmax = 80;
optionslam1.truesol = double(mp(lam1)\mp(b1));

% Run exact and double precision HSCG (2-term)
resultslam1 = pcge(lam1, b1, eye(size(lam1,1)), tol, optionslam1);
resultslam1d = pcgd(double(lam1), double(b1), eye(size(lam1,1)), tol, optionslam1);

% Run exact and double precision STCG (3-term)
[ xe,re,Aerre] = STCG_mp( lam1,b1,mp(zeros(n,1),68),nmax,68);
[ xd,rd,Aerrd] = STCG( lam1,b1,zeros(n,1),nmax);

% Cut-offs for plotting
i = find(resultslam1.error_A_norm < tol,1,'first');
if(~isempty(i))
    resultslam1.error_A_norm = resultslam1.error_A_norm(1:i);
end

i = find(Aerre < tol,1,'first');
if(~isempty(i))
    Aerre = Aerre(1:i);
end

i = find(resultslam1d.error_A_norm < tol,1,'first');
if(~isempty(i))
    resultslam1d.error_A_norm = resultslam1d.error_A_norm(1:i);
end

i = find(Aerrd < tol,1,'first');
if(~isempty(i))
    Aerrd = Aerrd(1:i);
end

% Plot convergence trajectories
figure()
semilogy(0:numel(resultslam1.error_A_norm)-1,resultslam1.error_A_norm,'r-','LineWidth',2);
hold on;
semilogy(0:numel(Aerre)-1,Aerre,'b--','LineWidth',2);
semilogy(0:numel(resultslam1d.error_A_norm)-1,resultslam1d.error_A_norm,'r:','LineWidth',2);
semilogy(0:numel(Aerrd)-1,Aerrd,'b:','LineWidth',2);

axis([0,40,1e-17,10])
set(gca,'FontSize',16)
print -depsc twovthree.eps