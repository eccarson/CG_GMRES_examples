%Generate diagonal "strakos matrix" test matrix but with cluster at high
%end
%n is dimension
%l1 and ln are smallest and largest eigenvalues, resp.
%0<p<=1 determines eigenvalue distribution
function A = strakosmatrix2(n, l1, ln, p)

d(1) = (l1);
d(n) = (ln);
for i = n-1:-1:2
    d(n-i+1) = ln - ((i-1)/(n-1))*(ln-l1)*(p)^(n-i);
end
A = diag(d);