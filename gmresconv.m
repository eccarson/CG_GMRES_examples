
n = 20;
mp.Digits(200);

% Problem with all eigenvalues = 1, but takes n iterations to converge

resnorms = mp([ones(1,n-1),1e-13,0]');
%B = mp(gallery('circul',[0,1,zeros(1,n-2)]));
evals = ones(n,1);%eig(B);

[A,b] = constructGMRESproblem(resnorms, evals);

[X,FLAG,RELRES,ITER,resvec] = gmres(mp(A),mp(b),n,1e-12,n);
%resvec = resvec./resvec(1);

figure()
semilogy(0:numel(resvec)-1,resvec, 'b-', 'LineWidth',2)
axis([0,20,1e-12,10])
set(gca,'FontSize',16)
print -depsc gmresconv1.eps


clear all
% Problem with cond number very high (+ evals equally spaced), but converges in 1 iteration
n = 20;

resnorms = mp([1,1e-13.*ones(1,n-1),0]');
evals = mp(diag(strakosmatrix(n, 1e-6, 1e10, 1)));

[A,b] = constructGMRESproblem(resnorms, evals);

[X,FLAG,RELRES,ITER,resvec] = gmres(mp(A),mp(b),n,1e-12,n);
%resvec = resvec./resvec(1);

figure()
semilogy(0:numel(resvec)-1,resvec, 'b-', 'LineWidth',2)
axis([0,20,1e-12,10])
set(gca,'FontSize',16)
print -depsc gmresconv2.eps