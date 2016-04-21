function [L,R] = Householder_implicit(A)

% function [L,R] = Householder_implicit(A)
%
% Computes the QR factorisation of A using Householder 
% transformations. The columns of the lower triangular
% matrix L will contain the vectors v of the subsequent
% Householder transformations.

[m,n] = size(A);

if m < n,
  disp('A cannot have more columns than rows.')
  return
end

R = A;
L = zeros(m,n);

for k = 1:n
   x = R(k:m,k);
   if x(1) == 0
       sgn = 1;
   else
       sgn = sign(x(1));
   end
   v = sgn*norm(x)*eye(m-k+1,1) + x;
   v = v/norm(v);
   L(k:m,k) = v;
   R(k:m,k:n) = R(k:m,k:n) - 2*v*(v'*R(k:m,k:n));
end

end
