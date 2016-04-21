function [Q,R] = Householder_explicit(A)

% [Q,R] = Householder_explicit(A)
%
% Computes the QR factorisation of A using Householder 
% transformations. Performs the orthogonal
% triangularisation explicitly.

[m,n] = size(A);

if m < n,
  disp('A cannot have more columns than rows.')
  return
end

% This will eventually contain the product Q = Q_1 Q_2 ... Q_n
Q = eye(m);

% This will eventually contain the orthogonal triangularisation of A:
% R = Q* A = Q_n ... Q_2 Q_1 A
R = A;

for k = 1:n
   % Calculate the Householder reflector matrix F.
   x = R(k:m,k);
   if x(1) == 0
       sgn = 1;
   else
       sgn = sign(x(1));
   end
   v = sgn*norm(x)*eye(m-k+1,1) + x;
   F = eye(m-k+1) - 2*(v*v')/(v'*v);
   
   % Construct the unitary matrix Q_k.
   Qk = zeros(m);
   Qk(1:k-1, 1:k-1) = eye(k-1);
   Qk(k:m, k:m) = F;
   
   % Calculate Q_1 Q_2 ... Q_k
   Q = Q*Qk;
   
   % Calculate Q_k ... Q_2 Q_1 A
   R = Qk*R;
end

end
