% This script generates all plots for Conjugate Gradient examples

% "Mathematical behavior of CG for different eigenvalue distributions" and
% "Computational behavior of CG for different eigenvalue distributions"
accumexp(0.6)
clusterlikefp(0.6)
spectraldensity
spectraldensityexact

% "Mathematical behavior of CG for problems with clustered eigenvalues"
clusterexp(0.6)
spectraldensityclustexact

% "Worst-case CG"
worst

% "Numerical behavior of CG on model problems"
wisehart
poissoncg

% "Convergence of preconditioned CG"
pcgvscg

% "Computational behavior of different CG implementations"
fp2v3

% % "Convergence of CG for different right-hand sides"
% cgRHS

% "Residual versus error and stopping criteria"
resvserror

% "CG with "average" initial residuals"
randx0