function loglike = rvbpcamv_loglikelihood(Y, W, CovW, X, CovX, mu, v_mu, ...
                                          a_tau, b_tau, a_U, b_U)
% loglike = rvbpcamv_loglikelihood(Y,W,CovW,X,CovX,mu,v_mu,a_tau,b_tau,a_U,b_U)
% loglike = rvbpcamv_loglikelihood(Y,Q)
%
% where Q is a struct with the posterior approximation distribution (Q)
% parameters as fields.
%
% The function returns the predictive density for non-nan values of Y given
% the (independent and approximate) posterior distributions. Some
% marginalizations (W,u,tau) are done using sampling methods.

if isstruct(W)
  if isfield(W, 'A')
    Q = alex2jaakko(W);
    Q.a_tau = 1e6;%Q.a_tau;
    Q.b_tau = Q.a_tau/W.V;%Q.b_tau;
  else
    Q = W;
  end
  W = Q.W;
  CovW = Q.CovW;
  X = Q.X;
  CovX = Q.CovX;
  mu = Q.mu;
  v_mu = Q.v_mu;
  a_tau = Q.a_tau;
  b_tau = Q.b_tau;
  if isfield(Q, 'a_U')
    a_U = Q.a_U;
    b_U = Q.b_U;
  else
    a_U = ones(size(Y));
    b_U = ones(size(Y));
  end
end

if numel(a_tau) == 1
  a_tau = repmat(a_tau, [rows(W), 1]);
end
if numel(b_tau) == 1
  b_tau = repmat(b_tau, [rows(W), 1]);
end

[m,n] = size(Y);
d = size(X,1);

N = 1:n;
nsamp = 10;

mean_y = zeros(nsamp,1);
var_y = zeros(nsamp,1);

% Sample W
w = zeros(m,d,nsamp);
for i=1:m
  w(i,:,:) = reshape(mnormrnd(W(i,:)',CovW(:,:,i),nsamp), [1,d,nsamp]);
end

% Sample mu
muh = normrnd(repmat(mu,1,nsamp), repmat(sqrt(v_mu),1,nsamp));
% Sample tau
tauh = gamrnd(repmat(a_tau(:),1,nsamp), repmat(1./b_tau(:),1,nsamp));

loglike = 0;
for j=1:n
  % Sample X
  x = mnormrnd(X(:,j),CovX(:,:,j),nsamp);
  for i=1:m

    if ~isnan(Y(i,j))
    
      % Sampled estimation of loglikelihood.
      for k=1:nsamp
        mean_y(k) = w(i,:,k)*x(:,k) + muh(i,k);
        var_y(k) = b_U(i,j)/a_U(i,j) *  1/tauh(i,k);
      end
      nu_y = 2*a_U(i,j);
    
      loglike = loglike + log(mean(studpdf(Y(i,j), mean_y, var_y, nu_y)));
    end
  end
end
    
% $$$ for i=1:m
% $$$   jmv = ~isnan(Y(i,:));
% $$$   % Sample tau
% $$$   tau = gamrnd(a_tau(i), 1/b_tau(i), nsamp, 1);
% $$$   % Sample W (because W*X is not Gaussian)
% $$$   w = mnormrnd(W(i,:)', CovW(:,:,i), nsamp);
% $$$   for j=N(jmv)
% $$$     
% $$$     % Sample u
% $$$     u = gamrnd(a_U(i,j), 1/b_U(i,j), nsamp, 1);
% $$$     
% $$$     % Sampled estimation of loglikelihood.
% $$$     mean_y = w'*X(:,j) + mu(i);
% $$$     %    mean_y = W(i,:)*X(:,j) + mu(i);
% $$$     for k=1:nsamp
% $$$       % TODO: IS THIS CORRECT???
% $$$       var_y(k) = w(:,k)'*CovX(:,:,j)*w(:,k) + v_mu(i) + 1 / (u(k)*tau(k));
% $$$     end
% $$$ %    loglike = loglike + log(mean(normpdf(Y(i,j), mean_y, sqrt(var_y))));
% $$$     
% $$$     % OR: sample W,X,tau and use Student-t distribution ?
% $$$     
% $$$     
% $$$     mean_y = w'*X(:,j) + mu(i);
% $$$     %    mean_y = W(i,:)*X(:,j) + mu(i);
% $$$     for k=1:nsamp
% $$$       % TODO: IS THIS CORRECT???
% $$$       var_y(k) = w(:,k)'*CovX(:,:,j)*w(:,k) + v_mu(i) + 1 / (u(k)*tau(k));
% $$$     end
% $$$     loglike = loglike + log(mean(mstud_lpdf(Y(i,j), mean_y, var_y, nu_y)));
% $$$     
% $$$               
% $$$   end
% $$$ end
