
clear all

mp.Digits(200);

%create matrix and rhs

n = 35;
A = strakosmatrix(n, 0.1, 100, 0.65);
b = ones(n, 1); b = b./norm(b);

xt = double(mp(A)\mp(b));
options.truesol = xt;
options.xlim = 2*n;


%double precision CG
% Run Lanczos/CG in double precision with NO reorthogonalization 
[Vd,T,bet] = tridiag_lan_double(A,b,options.xlim,0);

% Recover CG solutions, residuals, and errors from Lanczos quantities
[Xd,Rd,resd] = comp_solutions_lanczos(T,Vd,bet);
[Rdt,resdt] = comp_trueres(A,b,Xd);
Xdt = comp_error(A,xt,Xd);
Xdt = Xdt./Xdt(1);

%"exact" CG
[Ve,T,bet] = tridiag_lan_exact(A,b,options.xlim);

% Recover CG solutions, residuals, and errors from Lanczos quantities
[Xex,Rex,resex] = comp_solutions_lanczos_exact(T,Ve,bet);
[Rext,resext] = comp_trueres_exact(A,b,Xex);
Xext = comp_error_exact(A,xt,Xex);
Xext = Xext./Xext(1);

lm1=60;
lm2=26;
figure
semilogy(0:lm2-1, Xext(1:lm2), 'b-','LineWidth',2)
hold on
semilogy(0:lm1-1, Xdt(1:lm1), 'bo','LineWidth',2, 'MarkerSize', 8);
axis([0,lm1,1e-16,10])
set(gca,'FontSize',16)
print -depsc cgtraj_dve.eps

%Compute k sequence
for i = 1:size(Vd,2)-1
    rnk(i) = rank(Vd(:,1:i), 1e-1);
end
mx = max(rnk);
for i = 1:mx
    kk(i) = find(rnk==i, 1, 'last');
end


%compute differences
Xdtk = Xdt(kk);
Xextk = Xext(1:mx);


rat = Xdtk./Xextk;
onem = abs(1-(Xdtk./Xextk));

lm2 = 26;
figure()
semilogy(0:lm2-1,Xext(1:lm2), 'b-','LineWidth',2)
hold on
semilogy(0:lm2-1, Xdtk(1:lm2), 'bo','LineWidth',2, 'MarkerSize', 8);
semilogy(0:lm2-1, rat(1:lm2), 'r:','LineWidth',2);
semilogy(0:lm2-1, onem(1:lm2), 'r--','LineWidth',2);

axis([0,lm2,1e-16,10])
set(gca,'FontSize',16)
print -depsc cgtraj.eps
