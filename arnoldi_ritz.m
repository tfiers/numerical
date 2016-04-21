function [H, Q, ritz_values] = arnoldi_ritz(A, b, maxit)

% function [H, Q] = arnoldi_ritz(A, b, maxit)
%
% Arnoldi iteratie
% 
% invoer:
% A     - ijle matrix
% b     - startvector
% maxit - aantal iteraties
%
% uitvoer:
% H     - Hessenberg matrix
% Q     - orthogonale matrix
% ritz_values  - Matrix met in kolom i de Ritz-waarden van stap i.

ritz_values = zeros(maxit);

Q(:,1) = b/norm(b);
for n=1:maxit
  v = A*Q(:,n);
  for j = 1:n
    H(j,n) = Q(:,j)'*v;
    v = v - H(j,n)*Q(:,j);
  end
  H(n+1,n) = norm(v);
  if H(n+1,n) <= 0, 
      break;
  end
  Q(:,n+1) = v/H(n+1,n);
  Hn = H(1:end-1,:);
  %ritz_values(1:n, n) = sort(real(eig(Hn)), 1, 'descend');
  ritz_values(1:n, n) = real(eig(Hn));
end;

