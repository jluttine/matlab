
function test_kalman_filter()

% Dynamical model (noisy resonator)
D = 2;
w = 0.2;
A = [cos(w) sin(w)/w
     -w*sin(w) cos(w)];
% $$$ A = [1 0
% $$$      -0.5 0.5];
q = randn(D);
Q = q*q';

% Observation model
M = 4;
C = randn(M,2);
r = randn(M);
R = 1e2*r*r';

% Initial state
x0 = zeros(D,1);
Covx0 = eye(D);
x = mvnrnd(x0,Covx0)';
%Covx0 = Covx0 * 1e10;

% Simulate data
N = 100;
Y = zeros(M,N);
X = zeros(D,N);
for n=1:N
  x = mvnrnd(A*x, Q)';
  Y(:,n) = mvnrnd(C*x, R)';
  X(:,n) = x;
end

%tsplot(Y)
%tsplot(X)
X_true = X;

% Kalman filter (method 1)
[X, CovX, logp] = kalman_filter(x0, Covx0, Y, A, Q, C, R, 1);
logp

varx = reshape(CovX, [D*D, N]);
varx = varx(1:(D+1):end,:);
tserrorplot(X', 2*sqrt(varx'))
addtsplot(X_true, 'r');

% Kalman filter (method 2)
[X, CovX, logp] = kalman_filter(x0, Covx0, Y, A, Q, C, R, 2);
logp

varx = reshape(CovX, [D*D, N]);
varx = varx(1:(D+1):end,:);
tserrorplot(X', 2*sqrt(varx'))
addtsplot(X_true, 'r');

% RTS smoother
[X, CovX, x0, Covx0, CovX_X, entropy] = rts_smoother(x0, Covx0, X, CovX, A, Q);

varx = reshape(CovX, [D*D, N]);
varx = varx(1:(D+1):end,:);
tserrorplot(X', 2*sqrt(varx'))
addtsplot(X_true, 'r');

entropy
