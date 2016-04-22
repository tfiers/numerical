function [x, itx] = NMB_gmres(A, b)

% function [x, itx] = NMB_gmres(A, b)
%
% Solves Ax=b for x using the generalised minimal residuals method.
% The columns of 'itx' contain the approximations for x through the
% iteration steps.



% Each iteration step consists of three stages:
%
%  1) Step n of Arnoldi iteration to determine H and Q.
%
%  2) Solve least-squares problem ||rn|| = || H y - ||b||e1 ||    for y.
%                 (In ch. 11:     ||r || = || A x -      b  ||    for x).
%
%     ( Dimensions: (n-1) x n ; Hessenberg structure.)
%     With QR-factorisation:
%           2.1) Compute reduced QR factorisation H = QQ*RR
%           2.2) Compute the vector QQ' * ||b||e1
%           2.3) Solve the upper triangular system RR y = QQ'*||b||e1 for y.
%
%     (Optionally: save further work by using an updating process to get 
%      QR-factorisation of H from that of the previous H, using a Givens
%      rotation).
%
%  3) Calculate xn = Q y


[m,~] = size(A);
maxit = ceil(1.2*m);
%maxit = 150;
itx = zeros(m,maxit);

warning off all;
Q(:,1) = b/norm(b);
for n=1:maxit
    
    % 1) Arnoldi iteration
    v = A*Q(:,n);
    for j = 1:n
      H(j,n) = Q(:,j)'*v;
      v = v - H(j,n)*Q(:,j);
    end
    H(n+1,n) = norm(v);
    if H(n+1,n) <= 0, 
        break;
    end
    
    % 2) Least squares problem.
    [QQ,RR] = qr(H);
    bb = QQ' * norm(b)*eye(n+1,1);
    y = RR\bb;
    
    % 3) Calculate x
    x = Q*y;
    itx(:,n) = x;
    
    % (Finally, expand Q).
    Q(:,n+1) = v/H(n+1,n); 
end
warning on all;

end
