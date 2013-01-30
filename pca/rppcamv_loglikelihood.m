function loglike = rppcamv_loglikelihood(Y, W, X, invB, mu, tau, a_u, b_u)

if nargin == 2
  Q = W;
  W = Q.W;
  X = Q.X;
  mu = Q.mu;
  tau = Q.tau;
  a_u = Q.a_u;
  b_u = Q.b_u;
  invB = Q.invB;
end

[m,n] = size(Y);
d = size(W,2);

N = 1:n;
loglike = 0;

if isempty(X)
  X = zeros(d,n);
end
if isempty(invB)
  invB = repmat(eye(d), [1,1,n]);
end

for j=1:n
  imv = ~isnan(Y(:,j));
  if any(imv)
    mean_y = W(imv,:)*X(:,j) + mu(imv);
    Cov_y = b_u(j)/a_u(j) * (W(imv,:)*invB(:,:,j)*W(imv,:)'+eye(sum(imv))/tau);
    nu_y = 2*a_u(j);
    loglike = loglike + mstud_lpdf(Y(imv,j), mean_y, Cov_y, nu_y);
  end
end

