%Generate diagonal "strakos matrix" test matrix
%n is dimension
%l1 and ln are smallest and largest eigenvalues, resp.
%0<p<=1 determines eigenvalue distribution
function A = strakosmatrix(n, l1, ln, p)

d(1) = (l1);
d(n) = (ln);
for i = 2:n-1
    d(i) = l1 + ((i-1)/(n-1))*(ln-l1)*(p)^(n-i);
end
A = diag(d);