function [A,b] = constructGMRESproblem(resnorms, evals)

n = numel(evals);


for i = 1:n
    g(i) = sqrt(resnorms(i)^2 - resnorms(i+1)^2 );
end

V = eye(n);
b = g';

B = [b, V(:, 1:n-1)];

syms x
poly = (x-evals(1));
for i = 2:numel(evals)
    poly = poly*(x-evals(i));
end


cvec = coeffs(poly);
AB = diag(ones(n-1,1),-1);
AB(:,n) = cvec(2:end);





A = mp(double(B))*mp(double(AB))*mp(inv(mp(double(B))));


end