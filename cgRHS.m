% This script creates plots for the example:
% "CG convergence for different right-hand sides"
%
% Joerg Liesen, February 13, 2020
%
% Note: requires sdpt3 library: https://github.com/sqlp/sdpt3


N=100; ee=ones(N,1); x0=zeros(N,1); tol=1e-12; MAXIT=N;
A = spdiags([-ee,2*ee,-ee],-1:1,N,N);
[Q,D] = eig(full(A));
x1=Q*ee; b1=A*x1; [Anorm1] = my_cg(A,b1,x0,x1,tol,MAXIT);  
b2=Q*ee; x2=A\b2; [Anorm2] = my_cg(A,b2,x0,x2,tol,MAXIT);  

lam = eig(full(A)); kappa = max(lam)/min(lam);
factor = (sqrt(kappa)-1)/(sqrt(kappa)+1);
out(1)=1; kappa_bd(1)=2;
for k=1:length(Anorm1)-1
    [blk,Avec,C,b,X0,y0,Z0,objval,p] = cheby0(lam,k,1);
    out(k+1)=objval;
    kappa_bd(k+1)=2*factor^k;
end

f1 = figure;
figure(f1), clf
plot(real(lam),imag(lam),'.','MarkerSize',8)
axis([-0.1 4.1 -0.1 0.1])
yticks([0.0])
yticklabels({'\lambda(A)'})
ax = gca;
set(ax,'FontSize',16)
%print -depsc evals-1d-poisson.eps

f2 = figure;
figure(f2), clf, 
semilogy([0:length(Anorm1)-1]',Anorm1,'b--','LineWidth',2), hold on
semilogy([0:length(Anorm2)-1]',Anorm2,'r--','LineWidth',2)
semilogy([0:length(out)-1],out,'k:','LineWidth',2)
semilogy([0:length(kappa_bd)-1],kappa_bd,'k','LineWidth',2)
axis([0 length(Anorm1)-1 0.5*1e-3 2])
ax = gca;
set(ax,'FontSize',14)
print -depsc conv-1d-poisson.eps

f3 = figure;
figure(f3), clf
semilogy([1:N]',abs(Q'*b1),'b--','LineWidth',2), 
hold on, semilogy([1:N]',abs(Q'*b2),'r--','LineWidth',2)
axis([1 N 1e-3 1e+1])
ax = gca;
set(ax,'FontSize',14)
print -depsc coeffs-1d-poisson.eps