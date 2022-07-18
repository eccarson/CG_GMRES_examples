% This script generates convergence plot for the example:
% "Convergence of preconditioned CG"

% Set digits to simulate exact arithmetic
mp.Digits(100);

% Construct diagonal matrix 
A = strakosmatrix(40, .001, 100, .1);
A = mp(A);

% Construct "preconditioned" matrix
B = diag(linspace(10,100,40));
B = mp(B);

% Construct preconditioner P
P = mp(A)/mp(B);

% Construct right-hand side
b = [ones(size(A,1),1)];
b = b./norm(b);

% Set CG options
tol = 1e-16;
options.xlim=50;
options.truesol = double(mp(A)\mp(b));

% Run exact CG for the unpreconditioned system
resultsnp = pcge(A, b, eye(size(A,1)), tol, options);


% Run exact CG for the preconditioned system
results1 = pcge(mp(P)\mp(A), mp(P)\mp(b), eye(size(A,1)), tol, options);

% Cut-offs for plotting
i = find(resultsnp.error_A_norm < 7e-16,1,'first');
resultsnp.error_A_norm = resultsnp.error_A_norm(1:i);

i = find(results1.error_A_norm < 7e-16,1,'first');
results1.error_A_norm = results1.error_A_norm(1:i);


% Plot error norms for both systems
figure()
semilogy(0:numel(resultsnp.error_A_norm)-1,resultsnp.error_A_norm,'b--','LineWidth',2);
hold on;
semilogy(0:numel(results1.error_A_norm)-1,results1.error_A_norm,'r--','LineWidth',2);

axis([0,40,1e-18,10])
set(gca,'FontSize',15)