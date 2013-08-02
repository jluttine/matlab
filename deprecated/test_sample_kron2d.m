function test_sample_kron2d
warning('This function is deprecated')


% Length-scale=30 and PCG-iterations=20
% 500x500 = 2.7s
% 1000x1000 = 8.9s
% 2000x2000 = 27.7s
% 4000x4000 = 120.2 s


N1 = 500;
x1 = 1:N1;
K1 = gp_cov_pp(log(20), x1, x1);
%I1 = 3e-1*speye(N1);

N2 = 500;
x2 = 1:N2;
K2 = gp_cov_pp(log(20), x2, x2);
%I2 = speye(N2);

% $$$ subplot(2,1,1)
% $$$ imagesc(K1)
% $$$ subplot(2,1,2)
% $$$ imagesc(K2)
% $$$ return


%
% Simulate 2-D data
%

s = 1e-0;
[X2,X1] = meshgrid(x2,x1);
freq1 = 0.03*2*pi;
freq2 = 0.07*2*pi;
% Y_noiseless = sin(freq1*X1).*sin(freq1*X2) +
% sin(freq2*X1).*sin(freq2*X2);
Y_noiseless = chol(K1)' * randn(N1,N2) * chol(K2,'lower')';
Y = Y_noiseless + s*randn(N1,N2);

clf
subplot(2,1,1)
%shading interp
%contourf(X1,X2,Y)
imagesc(Y)


% $$$ % Use preconditioner kron(K1+s^2*I, K2+s2^I) ?
% $$$ C = kron(K1,K2)+s^2*speye(N1*N2);
% $$$ condest(C)
% $$$ C_precond = kron(K1+s*speye(N1), K2+s*speye(N2));
% $$$ condest(C_precond\C)
% $$$ return

%
% Draw samples from the posterior using C-G
%

M = 1; % number of samples

O = zeros(N1,N1);

% $$$ L1 = chol(K1, 'lower');
% $$$ U2 = chol(K2, 'upper');
%[LD1,p,q] = ldlchol(K1);
%nnz_LD1 = nnz(LD1)
LD1 = ldlchol(K1);
%nnz_LD1 = nnz(LD1)
[L,D] = ldlsplit(LD1);
L1 = L*sqrt(D);

LD2 = ldlchol(K2);
[L,D] = ldlsplit(LD2);
L2 = L*sqrt(D);


%afun = @(x) kronprod_for_cg(K1,K2,s^2,x);

C1 = K1+s*speye(N1);
C2 = K2+s*speye(N2);
LD_C1 = ldlchol(C1);
LD_C2 = ldlchol(C2);
%precondfun = [];
%precondfun = @(x) kronsolve_for_cg2(LD_C1, LD_C2, x);
% precondfun = @(x) kronsolve_for_cg1(C1, C2, x)

Kfun = @(x) (kronprod(K2,K1,reshape(x,N1,N2),'vector') + s^2*x);
K1fun = @(x) kronprod(K2,K1,reshape(x,N1,N2),'vector');
L1fun = @(x) kronprod(L2,L1,reshape(x,N1,N2),'vector');
L2fun = @(x) s*x;
Mfun = @(x) linsolve_kron_ldlchol(LD_C2, LD_C1, reshape(x,N1,N2), 'vector');

tic
for m=1:M
  F = gaussian_rand_pcg(Y(:), Kfun, K1fun, L1fun(randn(N1,N2)), L2fun(randn(N1*N2,1)), Mfun);
  F = reshape(F,N1,N2);
% $$$   Z0 = randn(N1,N2);
% $$$   Z = randn(N1,N2);
% $$$   % F0 = kronprod(L1,U2,Z0);
% $$$   F0 = kronprod(L1,L2,Z0);
% $$$   G = F0 + s*Z - Y;
% $$$   X = reshape( pcg(afun, G(:), 1e-4, 20, precondfun), N1, N2 );
% $$$   F = F0 - kronprod(K1,K2,X);
end
toc

subplot(2,1,2)
%shading interp
%contourf(X1,X2,F)
imagesc(F)

rmse_Y = rmse(Y_noiseless,Y)
rmse_F = rmse(Y_noiseless,F)


%
% Brute-force
%

% $$$ nnz_K1 = nnz(K1)
% $$$ nnz_K2 = nnz(K2)
% $$$ 
% $$$ K = kron(K2,K1);
% $$$ tic
% $$$ LD = ldlchol(K);
% $$$ toc


return

mu1 = mean(f1,2);
std1 = std(f1,1,2);

%
% Draw samples from the brute force posterior
%

tic
mu = K1*inv(C1)*y1;
V = K1 - K1*inv(C1)*K1;
f2 = bsxfun(@plus, mu, chol(V,'lower')*randn(N1,M));
toc

hold on
plot(x1,f2, 'r')

mu2 = mean(f2,2);
std2 = std(f2,1,2);

subplot(2,1,2)

plot(x1, [mu1, mu1-2*std1, mu1+2*std1], 'b');
hold on
plot(x1, [mu2, mu2-2*std2, mu2+2*std2], 'r');

% $$$ function x = kronprod_for_cg(K1,K2,s2,x);
% $$$ x = reshape(x,size(K1,1),size(K2,1));
% $$$ x = K1*x*K2 + s2.*x;
% $$$ x = x(:);
%afun = @(x) reshape(kronprod(K1,K2,reshape(x,N1,N2))+s^2*, N1*N2, 1);

% $$$ function Y = kronprod(A, B, X)
% $$$ % Evaluates: reshape( kron(B, A) * X(:), size(X) )
% $$$ Y = A*X*B';

% $$$ function x = kronsolve_for_cg2(LD1, LD2, x);
% $$$ x = reshape(x,size(LD1,1),size(LD2,1));
% $$$ x = ldlsolve(LD1,x);
% $$$ x = ldlsolve(LD2,x')';
% $$$ %x = (K1\x)/K2;
% $$$ x = x(:);
% $$$ 
% $$$ function x = kronsolve_for_cg1(K1, K2, x);
% $$$ x = reshape(x,size(K1,1),size(K2,1));
% $$$ x = (K1\x)/K2;
% $$$ x = x(:);
