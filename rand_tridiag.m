function A = rand_tridiag(n)

% function A = rand_tridiag(n)
%
% Generates a random tridiagonal matrix of dimension n.

A = diag(rand(n,1));
A = A + diag(rand(n-1,1),1);
A = A + triu(A,1)';

end
