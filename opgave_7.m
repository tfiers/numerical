
density = 1e-3;
%A = sprand(1000, 1000, density);
%b = rand(1000, 1);
maxit = 100;

[H,Q,ritz_values] = arnoldi_ritz(A, b, maxit);

pl = newplotlist;
x = 1:maxit;
for i = 1:maxit
    pl = addplotlist(pl, '', x, ritz_values(i, :), 'b.');
end

real_eigs = real(eigs(A,200));

%interesting = real_eigs(abs(real_eigs) > 1e-4);
interesting = real_eigs;

[num,~] = size(interesting);
x = maxit*ones(num,1);
pl = addplotlist(pl, '', x, interesting, 'ro');

doplotlist(pl, 'plot');
