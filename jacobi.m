function [V,D] = jacobi(A, tol)

% function [V,D] = jacobi(A, tol)
%
% Computes the eigenvalue factorisation V*D*V' of a symmetric matrix A,
% up to a given tolerance level, using the Jacobi method.

[n,m] = size(A);
if n~=m,
  disp('A is not square.')
  return
end
if n<2
  disp('A needs to have at least dimension 2.')
  return
end

MAXIT = 100;

V = eye(n);
for it = 1:MAXIT
    for i = 1:n-1
        for j = i+1:n
            theta = atan(2*A(i,j)/(A(j,j)-A(i,i)))/2;
            c = cos(theta);
            s = sin(theta);
            J = eye(n);
            J([i,j], [i,j]) = [c s; -s c];
            A = J'*A*J;
            V = V*J;
        end
    end
    if norm(A-diag(diag(A))) / max(diag(A)) < tol
        disp(it);
        break
    end
end
D = A;

end
