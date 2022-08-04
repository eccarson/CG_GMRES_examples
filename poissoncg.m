% This script generates convergence plots for the example:
% "Numerical behavior of CG on a Poisson model problem"

% Construct matrix and RHS
N = 50;
A = gallery('poisson', N);
P = speye(N^2);
b = ones(N^2,1);

% Set options for CG
options.xlim = 120;
tol = 1e-16;
options.truesol = mp(A,64)\mp(b,64);
xt = options.truesol;

%exact CG

% Run Lanczos/CG in double precision with reorthogonalization (simulates
% exact arithmetic)
[Ve,T,bet] = tridiag_lan_double(A,b,options.xlim,2);

% Compute loss of orthogonality in Krylov basis
for i = 1:size(Ve,2)
    looe(i) = norm(eye(i) - Ve(:,1:i)'*Ve(:,1:i), 'fro');
end

% Recover CG solutions, residuals, and errors from Lanczos quantities
[Xex,Rex,resex] = comp_solutions_lanczos(T,Ve,bet);
[Rext,resext] = comp_trueres(A,b,Xex);
Xext = comp_error(A,xt,Xex);
Xext = Xext./Xext(1);

% Run Lanczos/CG in double precision with NO reorthogonalization 
[Vd,T,bet] = tridiag_lan_double(A,b,options.xlim,0);

% Compute loss of orthogonality in Krylov basis
for i = 1:size(Vd,2)
    lood(i) = norm(eye(i) - Vd(:,1:i)'*Vd(:,1:i), 'fro');
end

% Recover CG solutions, residuals, and errors from Lanczos quantities
[Xd,Rd,resd] = comp_solutions_lanczos(T,Vd,bet);
[Rdt,resdt] = comp_trueres(A,b,Xd);
Xdt = comp_error(A,xt,Xd);
Xdt = Xdt./Xdt(1);

% Plot loss of orthogonality and A-norms of the errors for both
% reorthogonalized case and case with no reorthogonalization
figure()
semilogy(1:numel(looe), looe, ':k','LineWidth',2);
hold on
semilogy(1:numel(Xext), Xext, '-.k','LineWidth',2);
semilogy(1:numel(lood), lood, '--k','LineWidth',2);
semilogy(1:numel(Xdt), Xdt, '-k','LineWidth',2);
axis([0,120,1e-16,1])
set(gca,'FontSize',16)
print -depsc poissonCG.eps