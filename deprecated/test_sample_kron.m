function test_sample_kron
warning('This function is deprecated')


N1 = 200;
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
% Draw samples from the posterior
%

M = 1000; % number of samples

z = randn(2*N1,M);
O = zeros(N1,N1);
L1 = chol(K1, 'lower');
f1 = [L1 O]*z - K1 * (C1 \ bsxfun(@minus, [L1 I1.^0.5]*z, y1));

hold on
plot(x1, f1, 'b')

mu1 = mean(f1,2);
std1 = std(f1,1,2);

%
% Draw samples from the brute force posterior
%

mu = K1*inv(C1)*y1;
V = K1 - K1*inv(C1)*K1;
f2 = bsxfun(@plus, mu, chol(V,'lower')*randn(N1,M));

hold on
plot(x1,f2, 'r')

mu2 = mean(f2,2);
std2 = std(f2,1,2);

subplot(2,1,2)

plot(x1, [mu1, mu1-2*std1, mu1+2*std1], 'b');
hold on
plot(x1, [mu2, mu2-2*std2, mu2+2*std2], 'r');