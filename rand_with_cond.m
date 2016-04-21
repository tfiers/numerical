function M = rand_with_cond(n, kappa)

% function M = rand_with_cond(n, kappa)
%
% Returns a random n by n matrix with the 
% specified condition kappa.

[Q,~] = qr(rand(n));
D = diag(ones(n,1)); D(1,1) = kappa;
[V,~] = qr(rand(n));

M = Q*D*V';

end

