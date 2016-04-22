maxit = 1000;

%n = 1000;
%A = rand(n);
n = 35;
A = mat1; % load mat1.txt

[V,D] = eig(A);
startvec = V(:,1); %rand(n,1)
startvec = startvec+0.1*rand(n,1);

%es_qr = qr_shiftrayleigh(A,maxit);
[~, es_qr] = qr_shiftrayleigh(A,maxit);
es_gi = gelijktijdige_it(A,rand(n),maxit);
approxs = rayleigh(A,startvec,maxit);
disp(approxs(end));
es_rqi = abs( approxs - max(eig(A)) );

x = 1:maxit;

%semilogy(x,es_rqi, x,es_gi, 1:size(es_qr,2),es_qr);
loglog(x,es_rqi, x,es_gi, 1:size(es_qr,2),es_qr);
%semilogy(x,es_gi, x,es_qr);
xlabel 'Iteratiestap'
ylabel 'Fout op grootste eigenwaarde'