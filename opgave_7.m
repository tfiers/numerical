
density = 2e-3;
A = sprand(1000, 1000, density);
b = rand(1000, 1);
maxit = 100;

[H,Q,ritz_values] = arnoldi_ritz(A, b, maxit);

pl = newplotlist;
x = 1:maxit;
for i = 1:maxit
    pl = addplotlist(pl, '', x, ritz_values(i, :), 'b.');
end

real_eigs = sort(real(eigs(A,200)), 'descend');

%interesting = real_eigs(abs(real_eigs) > 1e-4);
interesting = real_eigs;

[num,~] = size(interesting);
x = maxit*ones(num,1);
pl = addplotlist(pl, '', x, interesting, 'ro');

close all;
figure;
doplotlist(pl, 'plot');
xlabel 'Iteratiestap';
ylabel 'Reëel deel Ritz- en eigenwaarden';

rrr = sort(ritz_values, 'descend');

% figure;
% semilogy(abs(rrr(1,:)-real_eigs(1)));
% xlabel 'Iteratiestap';
% ylabel 'Fout op de grootste eigenwaarde';
% 
% figure;
% semilogy(abs(rrr(end,:)-real_eigs(end)));
% xlabel 'Iteratiestap';
% ylabel 'Fout op de kleinste eigenwaarde';

figure;
x = 1:maxit;
e1 = abs(rrr(1,:)-real_eigs(1));
e2 = abs(rrr(2,:)-real_eigs(2));
e3 = abs(rrr(end-1,:)-real_eigs(end-1));
e4 = abs(rrr(end,:)-real_eigs(end));
semilogy(x,e1, x,e2, x,e3, x,e4);
xlabel 'Iteratiestap';
ylabel 'Fout op eigenwaarden';

