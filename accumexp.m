% This function, run with rho=0.6, generates the convergence plots for the examples:
% 1) "Mathematical behavior of CG for different eigenvalue distributions"
% 2) "Computational behavior of CG for different eigenvalue distributions"

function accumexp(rho)

n = 30;
l1 = .1;
ln = 1e3;

% Construct 3 different eigenvalue distributions
lam1 = diag(strakosmatrix(n, l1, ln, rho));
lam2 = diag(strakosmatrix2(n, l1, ln, rho));
lam3 = diag(strakosmatrix(n, l1, ln, 1));


% Plot eigenvalue distributions
figure()
plot(real(lam1),imag(lam1)-0.0,'r.','MarkerSize',16), hold on
plot(real(lam2),imag(lam2)+0.5,'b.','MarkerSize',16), hold on
plot(real(lam3),imag(lam3)-0.5,'g.','MarkerSize',16), hold on
axis([-100,1100,-1,1])
%yticks([-.5,0,.5])
%yticklabels({'\Lambda(A_3)','\Lambda(A_2)','\Lambda(A_1)'})
set(gca,'YTickLabel',[]);
set(gca, 'fontsize', 16);
print -depsc eigdist.eps

% Construct diagonal matrices and right-hand sides
lam1 = diag(lam1); lam2 = diag(lam2); lam3 = diag(lam3);
b = ones(size(lam1,1),1);
b = b./norm(b);

% CG tolerance
tol = 1e-14;

% axis limit parameters for plotting
optionslam1.xlim=80; 
optionslam2.xlim=80;
optionslam3.xlim=60;

% Run exact CG for the 3 systems
optionslam1.truesol = lam1\b;
resultslam1 = pcge(lam1, b, eye(size(lam1,1)), tol, optionslam1);

optionslam2.truesol = lam2\b;
resultslam2 = pcge(lam2, b, eye(size(lam2,1)), tol, optionslam2);

optionslam3.truesol = lam3\b;
resultslam3 = pcge(lam3, b, eye(size(lam3,1)), tol, optionslam3);

% Run double precision CG for the 3 systems
resultslam1d = pcgd(double(lam1), double(b), eye(size(lam1,1)), tol, optionslam1);
resultslam2d = pcgd(double(lam2), double(b), eye(size(lam2,1)), tol, optionslam2);
resultslam3d = pcgd(double(lam3), double(b), eye(size(lam3,1)), tol, optionslam3);

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

i = find(resultslam1d.error_A_norm < tol,1,'first');
if(~isempty(i))
    resultslam1d.error_A_norm = resultslam1d.error_A_norm(1:i);
end

i = find(resultslam2d.error_A_norm < tol,1,'first');
if(~isempty(i))
    resultslam2d.error_A_norm = resultslam2d.error_A_norm(1:i);
end

i = find(resultslam3d.error_A_norm < 7e-16,1,'first');
if(~isempty(i))
    resultslam3d.error_A_norm = resultslam3d.error_A_norm(1:i);
end

% Plot convergence for exact CG
figure()
semilogy(0:numel(resultslam2.error_A_norm)-1,resultslam2.error_A_norm,'b-','LineWidth',2);
hold on;
semilogy(0:numel(resultslam1.error_A_norm)-1,resultslam1.error_A_norm,'r-','LineWidth',2);
semilogy(0:numel(resultslam3.error_A_norm)-1,resultslam3.error_A_norm,'g-','LineWidth',2);

axis([0,50,1e-8,10])
%legend( 'acc. to the right', 'acc. to the left', 'equally spaced','Interpreter','latex');
set(gca,'FontSize',16)
print -depsc exactcgdist.eps

% Plot convergence for double precision CG
figure()
semilogy(0:numel(resultslam2d.error_A_norm)-1,resultslam2d.error_A_norm,'b-','LineWidth',2);
hold on;
semilogy(0:numel(resultslam1d.error_A_norm)-1,resultslam1d.error_A_norm,'r-','LineWidth',2);
semilogy(0:numel(resultslam3d.error_A_norm)-1,resultslam3d.error_A_norm,'g-','LineWidth',2);

axis([0,50,1e-8,10])
%legend('acc. to the right', 'acc. to the left', 'equally spaced','Interpreter','latex');
set(gca,'FontSize',16)
print -depsc dpcgdist.eps
%saveas(gcf,'fpcg_new.pdf')



