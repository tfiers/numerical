function A = rand_sym(n)

% function A = rand_sym(n)
%
% Generates a random symmetric matrix of dimension n.

A = triu(rand(n));
A = A + triu(A,1)';

end

