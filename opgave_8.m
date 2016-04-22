m = 100;
alpha = 100; % 1, 5, 10, 100
A = sprand(m,m,0.5);
A = A + alpha*speye(m); A=A/norm(A,1);
b = rand(m,1);

[x,itx] = NMB_gmres(A,b);
maxit = size(itx,2);

sol = A\b;

res_norms = zeros(1,maxit);
errors = zeros(1,maxit);
for i = 1:maxit
    res_norms(i) = norm(A*itx(:,i) - b);
    errors(i) = norm(itx(:,i) - sol);
end

ix = 1:maxit;

close all;
figure;
plot(ix,res_norms, ix,errors);
xlabel 'Iteratiestap'

figure;
semilogy(ix,res_norms, ix,errors);
xlabel 'Iteratiestap'

