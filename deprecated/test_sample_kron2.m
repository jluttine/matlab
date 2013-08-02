function test_sample_kron2
warning('This function is deprecated')


N1 = 1000;
x1 = 1:N1;
K1 = gp_cov_pp(log(30), x1, x1);
I1 = 3e-1*speye(N1);

N2 = 100;
x2 = 1:N2;
K2 = gp_cov_pp(log(10), x2, x2);
I2 = speye(N2);

%
% Simulate 1-D data
%

C1 = K1 + I1;
L_C1 = chol(C1, 'lower');
y1 = L_C1 * randn(N1,1);

clf
subplot(2,1,1)

%figure
plot(x1,y1,'r+')

%
% Draw samples from the posterior using C-G
%

M = 100; % number of samples

% afun = @(x) reshape(kronprod(K1+I1,K2+I2,reshape(x,N1,N2)), N1*N2, 1);

O = zeros(N1,N1);
L1 = chol(K1, 'lower');
f1 = zeros(N1,M);
tic
for m=1:M
  z = randn(2*N1,1);
  f1(:,m) = [L1 O]*z - K1 * pcg(C1, [L1 I1.^0.5]*z - y1);
end
toc

hold on
plot(x1, f1, 'b')

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