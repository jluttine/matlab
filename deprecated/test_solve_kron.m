
function test_solve_kron
warning('This function is deprecated')


N1 = 150;
x1 = 1:N1;
K1 = gp_cov_pp(log(10), x1, x1);
I1 = speye(N1);

N2 = 100;
x2 = 1:N2;
K2 = gp_cov_pp(log(10), x2, x2);
I2 = speye(N2);

% Full matrix approach
%K = kron(K1, K2);
%K_noise = kron(K1+I1, K2+I2);
%K_posterior = K - K * (K_noise \ K);

%
% Solve: kron(K1+I1,K2+I2) \ kron(K1,K2)
%

%K = kron(K2,K1);
%X_full = K * (kron(K2+I2,K1+I1) \ K);

% Kronecker approach
tic
C1 = K1 * ((K1 + I1) \ K1);
C2 = K2 * ((K2 + I2) \ K2);
X_kron = diag(C1) * diag(C2)';
%X_kron = dot(K1,C1,1)' * dot(K2,C2,1);
toc

% Conjugate gradient approach
%
% Computing the variances is veeery expensive...
tic
afun = @(x) reshape(kronprod(K1+I1,K2+I2,reshape(x,N1,N2)), N1*N2, 1);
X_cg = zeros(N1,N2);
for n1=1:N1
  for n2=1:N2
    fprintf('Progress: %d %%\n', floor(100*n1/N1));
    x = pcg(afun, kron(K2(:,n2),K1(:,n1)), [], 2);
    X_cg(n1,n2) = K1(n1,:) * reshape(x,N1,N2) * K2(:,n2);
  end
end
toc

full([X_kron(:), X_cg(:)])
%full([X_kron(:), X_cg(:), diag(X_full)])

function Y = kronprod(A, B, X)
% kron(B', A)
Y = A*X*B;