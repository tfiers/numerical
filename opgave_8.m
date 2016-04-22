res_norms = zeros(4,maxit);
errors = zeros(4,maxit);

m = 100;


alphas = [1, 5, 10, 100];

for a = 1:4
    alpha = alphas(a); % 1, 5, 10, 100
    
    A = sprand(m,m,0.5);
    b = rand(m,1);
    A = A + alpha*speye(m); A=A/norm(A,1);

    [x,itx] = NMB_gmres(A,b);
    maxit = size(itx,2);

    sol = A\b;

    for i = 1:maxit
        res_norms(a,i) = norm(A*itx(:,i) - b);
        errors(a,i) = norm(itx(:,i) - sol);
    end
end

ix = 1:maxit;

close all;
figure;
semilogy(ix,res_norms(4,:), 'b-', ...
         ix,   errors(4,:), 'b--',...
         ix,res_norms(3,:), 'g-', ...
         ix,   errors(3,:), 'g--',...
         ix,res_norms(2,:), 'r-', ...
         ix,   errors(2,:), 'r--',...
         ix,res_norms(1,:), 'c-', ...
         ix,   errors(1,:), 'c--');
xlabel 'Iteratiestap'

