function y = Apply_Q(L,b)

% y = Apply_Q(L,b)
%
% Computes y = Q*b, with Q from the QR factorisation
% of a matrix A. L is the lower triangular matrix 
% obtained from the 'Householder_implicit' factorisation 
% function applied to A.

[m,n] = size(L);

for k = 1:n
    v = L(k:m,k);
    b(k:m) = b(k:m) - 2*v*(v'*b(k:m));
end

y = b;

end
