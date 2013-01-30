
function test_rotationcost_gaussian

% WHY IS THE FIRST ELEMENT OF THE GRADIENT -10 ??? Add 1e-6 to R(1) and the
% gradient becomes zero.. Ah, fixed it by setting appropriate initialization
% for the posterior covariance (Cov(1)=1e-6) :) No problems thus! :)
% Hmm.. very small prior variance seems to be still a problem..?

D = 5;

% Fixed prior
Cov0 = diag([1e-6; ones(D-1,1)]);;
mu0 = [1; zeros(D-1,1)];

% Random posterior initialization
N = 10;
z = randn(D);
Cov = z*z';
Cov(1,:) = 0;
Cov(:,1) = 0;
Cov(1,1) = Cov0(1);
X = randn(D,N);
X(1,:) = 1;
XX = X*X' + N*Cov;


% Consider transformation mu <- R*mu and Cov <- R*Cov*R'
R = eye(D);

% $$$ c = [];
% $$$ for r1=20:-1:7
% $$$   R(1,1) = 1 - 10^(-r1);
% $$$   c = [c cost(R(:))];
% $$$   R(1,1) = 1 + 10^(-r1);
% $$$   c = [cost(R(:)) c];
% $$$ end
% $$$ figure
% $$$ plot(c)
% $$$ return

R(1) = 1;
cost(R(:));
%return

%R(1,1) = 1-eps;
mycheckgrad(@cost, R(:), 1e-6);
r = minimize(R(:), @cost, 1000);
R = reshape(r,[D,D])

% Update variables
X = R*X
XX = R*XX*R'
mean_X = mean(X,2)
CovX = 1/N*XX - mean_X*mean_X'

  function [f,df] = cost(r)
  
  R = reshape(r,[D,D]);
  
  logdet_Cov = 2*N*logdet(R);
  x_invCov0_x = traceprod(R*XX*R',inv(Cov0));
  x_invCov0_mu = sum(X,2)'*R'*inv(Cov0)*mu0;
  mu_invCov0_mu = 0;
  logdet_Cov0 = 0;
  grad_logdet_Cov = 2*N*inv(R)';
  grad_x_invCov0_x = 2*inv(Cov0)*R*XX;
  grad_x_invCov0_mu = inv(Cov0)*mu0*sum(X,2)';
  grad_mu_invCov0_mu = 0;
  grad_logdet_Cov0 = 0;
  
  [KL,dKL] = vb_rotationcost_gaussian(logdet_Cov, ...
                                      x_invCov0_x, ...
                                      x_invCov0_mu, ...
                                      mu_invCov0_mu, ...
                                      logdet_Cov0, ...
                                      grad_logdet_Cov, ...
                                      grad_x_invCov0_x, ...
                                      grad_x_invCov0_mu, ...
                                      grad_mu_invCov0_mu, ...
                                      grad_logdet_Cov0);
  
  f = KL;
  df = dKL(:);
  end
  
  
end
