function [V,D,errors] = jacobi(A, tol)

% function [V,D] = jacobi(A, tol)
%
% Computes the eigenvalue factorisation V*D*V^T of a symmetric matrix A,
% up to a given tolerance level.

[n,m] = size(A);
if n~=m,
  disp('A is not square.')
  return
end
if n<2
  disp('A needs to have at least dimension 2.')
  return
end

MAXIT = 10;
errors = [];
eigs = sort(eig(A));

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
%     disp(['Iteratiestap ' num2str(it)]);
%     disp(['Grootste diagonaalelement:      ' num2str(max(diag(A)))]);
%     disp(['Grootste niet-diagonaalelement: ' num2str(max(max(A-diag(diag(A)))))]);
%     disp(['Norm A met diagonaal 0:         ' num2str(norm(A-diag(diag(A))))]);
%     disp(' ');
    errors = [errors norm( sort(diag(A)) - eigs )];
    if norm(A-diag(diag(A))) / max(diag(A)) < tol
        %break
    end
end
D = A;

end
