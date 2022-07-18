% This script generates all plots for Conjugate Gradient examples

% "Mathematical behavior of CG for different eigenvalue distributions" and
% "Computational behavior of CG for different eigenvalue distributions"
accumexp(0.6)
spectraldensity
spectraldensityexact

% "Worst-case CG"
worst
wisehart

% "Numerical behavior of CG on a Poisson model problem"
poissoncg

% "Convergence of preconditioned CG"
pcgvscg

% "Computational behavior of different CG implementations"
fp2v3

% "Convergence of CG for different right-hand sides"
cgRHS

% "Residual versus error and stopping criteria"
resvserror

% "CG with "average" initial residuals"
randx0