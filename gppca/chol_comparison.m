%function chol_comparison
% Test what is the best way to evaluate Cholesky(-like) decomposition.

n = 1000;

% Large and ill-conditioned covariance matrix
C = gpcov(log(100), 1:n, 1:n);

% Regular Cholesky with diagonal addition
tic
Lchol = chol(C + 100*eps*eye(n), 'lower');
Lchol(1:10, 1:10)
toc

% LDL with negative values to zero
tic
[Lldl, Dldl] = ldl(C);
%Dldl(Dldl<0) = 0;
%Lldl = Lldl * sqrt(Dldl);
d = diag(Dldl);
d(d<0) = 0;
Lldl(1:10,1:10)
Lldl = bsxfun(@times, Lldl, sqrt(d(:))');
toc

% My LDL with negative values to zero
tic
[Lmyldl, Dmyldl] = ldl(C);
%Dldl(Dldl<0) = 0;
%Lldl = Lldl * sqrt(Dldl);
d = diag(Dmyldl);
d(d<0) = 0;
Lmyldl(1:10,1:10)
Lmyldl = bsxfun(@times, Lmyldl, sqrt(d(:))');
toc
%return

% Matlab's cholcov based on eigenvalue decomposition
% Note that this solution does not, in general, give triangular matrix!!!
% Also, this is very slow!
tic
Lcholcov = cholcov(C)';
Lcholcov = [Lcholcov, zeros(n, n-cols(Lcholcov))];
toc


% Errors:
error_chol = norm(C - Lchol*Lchol')
error_ldl = norm(C - Lldl*Lldl')
error_myldl = norm(C - Lmyldl*Lmyldl')
%error_ldl = norm(C - Lldl*Dldl*Lldl')
error_cholcov = norm(C - Lcholcov*Lcholcov')

% What other measures could be used?
% How accurate these methods are in solving linear equations?

x = randn(n,1);
y = C * x;

opts.LT = true;
opts.TRANSA = false;
x_chol = linsolve(Lchol, y, opts);
opts.TRANSA = true;
x_chol = linsolve(Lchol, x_chol, opts);

opts.LT = true;
opts.TRANSA = false;
x_ldl = linsolve(Lldl, y, opts);
opts.TRANSA = true;
x_ldl = linsolve(Lldl, x_ldl, opts);

opts.LT = true;
opts.TRANSA = false;
x_myldl = linsolve(Lmyldl, y, opts);
opts.TRANSA = true;
x_myldl = linsolve(Lmyldl, x_myldl, opts);

opts.LT = false;
opts.TRANSA = false;
x_cholcov = linsolve(Lcholcov, y, opts);
opts.TRANSA = true;
x_cholcov = linsolve(Lcholcov, x_cholcov, opts);

error_chol = norm(x - x_chol)
error_ldl = norm(x - x_ldl)
error_myldl = norm(x - x_myldl)
error_cholcov = norm(x - x_cholcov)
