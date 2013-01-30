
function test_mcmc_reflective()

N = 2;
mu = zeros(N,1);
Q = [1 1; -1 1]; %orth(randn(N));
D = (1-0.95)*eye(N);
D(1) = 1;
Cov = Q*D*Q';
%Cov = W'*W;
L = chol(Cov, 'lower');


get_logpdf = @(x) (@()gaussian_logpdf(x, mu, L));
get_dlogpdf = @(x) (@()gaussian_dlogpdf(x, mu, L));

x_init = 10*ones(N,1);
%x_init = zeros(N,1);

M = 1000;

% Test metropolis-hastings
y_mh = zeros(N,M);
if 1
  q = @(x) normrnd(x,0.3);
  logpdf_q = @(x,x0) normal_logpdf(x,x0,0.3);
  sampler = mcmc_init_metropolishastings(x_init, get_logpdf, q, logpdf_q);
  for m=1:M
    y_mh(:,m) = sampler();
  end
end

% Test hamiltonian
y_h = zeros(N,M);
if 1
  sampler = mcmc_init_hamiltonian(x_init, get_logpdf, get_dlogpdf, 0.3, 10);
  for m=1:M
    y_h(:,m) = sampler();
  end
end


% Test slice
y_s = zeros(N,M);
if 1
  sampler = mcmc_init_slicesampling(x_init, get_logpdf);
  for m=1:M
    y_s(:,m) = sampler();
  end
end

% Test inside reflective
y_r = zeros(N,M);
if 1
  sampler = mcmc_init_reflective(x_init, get_logpdf, get_dlogpdf, 0.01, ...
                                 300, 'type', 'inside');
  for m=1:M
    y_r(:,m) = sampler();
  end
  % TODO: Is this correct?
  y_r(:,isnan(y_r(1,:))) = [];
end

y = bsxfun(@plus, mu, L*randn(N,M));

figure(1)
clf();
plot(y(1,:), y(2,:), 'k.')
hold on
plot(y_mh(1,:), y_mh(2,:), 'r.')
plot(y_h(1,:), y_h(2,:), 'c.')
plot(y_s(1,:), y_s(2,:), 'b.')
plot(y_r(1,:), y_r(2,:), 'g.')

% Burn-in
M0 = 100;
y_mh = y_mh(:,M0:end);
y_h = y_h(:,M0:end);
y_s = y_s(:,M0:end);
y_r = y_r(:,M0:end);

M0 = 100;

T = Cov + mu*mu'
YY = y*y' / (size(y,2))
YY_mh = y_mh*y_mh' / (size(y_mh,2))
YY_h = y_h*y_h' / (size(y_h,2))
YY_s = y_s*y_s' / (size(y_s,2))
YY_r = y_r*y_r' / (size(y_r,2))

% Auto-correlation
figure(2)
clf();
lag = 100;
subplot(5,1,1)
plot(0:lag, acorr(y',lag))
subplot(5,1,2)
plot(0:lag, acorr(y_mh',lag))
subplot(5,1,3)
plot(0:lag, acorr(y_h',lag))
subplot(5,1,4)
plot(0:lag, acorr(y_s',lag))
subplot(5,1,5)
plot(0:lag, acorr(y_r',lag))

