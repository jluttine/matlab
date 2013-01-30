%
% [f, Covf, pseudoX, logtheta, Kpp, Kxp, Lp] = 
%     gplearn(y, V, X, pseudoX, kX, logtheta, ...)
%
% V is inverse covariance (i.e. precision) matrix for observation noise -
% sparse format is recommended.
%
% In practice, inv(V) should be diagonal. However, inv(V) can be more
% complex sparse matrix as long as it has zeros where kX is non-zero. That
% is, observation noise can be correlated if the corresponding latent
% function values are independent a priori. This more general property is
% used in PCA model, where inv(V) and kX are block-diagonal (although rows
% and columns may be in mixed order!) but the blocks "overlap" only on the
% diagonal, that is, they both have non-zero elements at the same time only
% on the diagonal.
%
% y = inv(V)*c, that is, c = V*y
% N( inv(V)*c | 0, inv(V) + Kxp inv(Kpp) Kpx )

function [f, Covf, pseudoX, logtheta, Kpp, Kxp, Lp] = gplearn(y, V, ...
                                                  X, pseudoX, kX, logtheta, ...
                                                  varargin)


opts = struct( ...
    'maxsearch', 10, ...
    'vy', [], ...
    'cholv', []);

[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end

m = cols(pseudoX); % number of pseudoinputs
d = rows(pseudoX); % dimensionality of the input space

if isempty(opts.cholv)
  LV = chol(V, 'lower');
else
  LV = opts.cholv;
end

if isempty(opts.vy)
  c = V*y;
else
  c = opts.vy;
end

if opts.maxsearch ~= 0
  % Learn pseudo-inputs and hyperparameters
  x0 = [pseudoX(:); logtheta(:)];
  func = @(x) bound2cost(x, c, V, LV, X, kX, m, d);
  [x, fiter] = minimize(x0, func, opts.maxsearch);
  pseudoX = reshape(x(1:(m*d)), [d,m]);
  logtheta = x((m*d+1):end);
end

% GP-prior covariances
if nargout >= 3
  Kpp = feval(kX, pseudoX, pseudoX, logtheta);
  Kpp = Kpp + 1e-6*eye(m);
end
if nargout >= 4
  Kxp = feval(kX, X, pseudoX, logtheta);
end
if nargout >= 5
  Lp = chol(Kpp, 'lower');
end

% Helpful variables
Lp = chol(Kpp, 'lower');
R = solve_tril(Lp, Kxp' * LV);
S = eye(m) + R*R';
LS = chol(S, 'lower');
z = LV \ c;
beta = solve_tril(LS, R*z);
invLsLp = solve_tril(LS, Lp');

% Posterior mean and covariance
f = invLsLp' * beta;
Covf = invLsLp' * invLsLp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, df] = bound2cost(x, c, V, LV, inputs, kX, m, d)

pseudoinputs = reshape(x(1:(m*d)), [d,m]);
logtheta = x((m*d+1):end);

[bound, dlogtheta, dpseudoinputs] = pseudobound(c, V, LV, inputs, pseudoinputs, logtheta, kX);

f = -bound;
df = -[dpseudoinputs(:); dlogtheta(:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bound, dlogtheta, dpseudoinputs] = ...
    pseudobound(c, V, LV, inputs, pseudoinputs, logtheta, kX)

d = rows(pseudoinputs); % dimensionality of pseudo-inputs
m = cols(pseudoinputs); % number of pseudo-inputs

% GP-prior (co)variances and their gradients
[Kpp, dKpp_logtheta, dKpp_pseudo] = feval(kX, pseudoinputs, pseudoinputs, logtheta);
[Kxp, dKxp_logtheta, dKxp_pseudo] = feval(kX, inputs, pseudoinputs, logtheta);
[Kx, dKx_logtheta] = feval(kX, inputs, [], logtheta);
Kpp = Kpp + 1e-6*eye(m);

% Helpful variables
Lp = chol(Kpp, 'lower');
R = solve_tril(Lp, Kxp' * LV);
S = eye(m) + R*R';
LS = chol(S, 'lower');
LA = Lp * LS;
z = LV \ c;
H = solve_tril(LS, R);
T = solve_triu(LA', H);
beta = H * z;
b = T * z;
invLpR = solve_triu(Lp', R);
diagV = full(diag(V));
dif = c - V*(Kxp*b);

% Lower bound for loglikelihood
bound = -logdettri(LS) - 0.5*z'*z + 0.5*beta'*beta ...
        - 0.5*Kx'*diagV ...
        + 0.5*traceprod(R,R');

% General gradient coefficients
Tpp = 0.5 * solve_triu(LA', solve_triu(LA', LS'*LS-eye(m))') ...
      - 0.5 * invLpR*invLpR';
Txp = -T * LV' ...
      + invLpR*LV';

% Gradient for hyperparameters
n = length(logtheta);
dlogtheta = zeros(n,1);
for i=1:n
  dlogtheta(i) = traceprod(dKpp_logtheta(:,:,i), Tpp) ...
      - 0.5 * b' * dKpp_logtheta(:,:,i) * b ...
      + traceprod(dKxp_logtheta(:,:,i), Txp) ...
      + dif' * dKxp_logtheta(:,:,i) * b ...
      - 0.5 * dKx_logtheta(:,:,i)' * diagV;
end

% $$$ size(Txp)
% $$$ size(dKxp_pseudo)

% Gradient for pseudo-inputs
dpseudoinputs = zeros(d,m);
Dpp = spalloc(m,m,2*m);
for i=1:m
  for j=1:d
    Dpp(:) = 0;
    Dpp(i,:) = dKpp_pseudo(j,:,i);
    Dpp(:,i) = dKpp_pseudo(j,:,i);
% $$$     Dxp = zeros(cols(inputs), m);
% $$$     Dxp(:,i) = dKxp_pseudo(j,:,i);

    dpseudoinputs(j,i) = traceprod(Dpp, Tpp) ...
        - 0.5 * b' * Dpp * b ...
        + dKxp_pseudo(j,:,i) * Txp(i,:)' ...
        + dif' * dKxp_pseudo(j,:,i)' * b(i);

% $$$     dpseudoinputs(j,i) = traceprod(Dpp, Tpp) ...
% $$$         - 0.5 * b' * Dpp * b ...
% $$$         + dKxp_pseudo(j,:,i) * Txp(i,:)' ...%traceprod(Dxp, Txp) ...
% $$$         + (c - V*Kxp*b)' * dKxp_pseudo(j,:,i)' * b(i);%(c - V*Kxp*b)' * Dxp * b;
% $$$ 
% $$$     dpseudoinputs(j,i) = ...%traceprod(dKpp_pseudo(:,:,i), Tpp) ...
% $$$         ...%- 0.5 * b' * dKpp_pseudo(:,:,i) * b ...
% $$$         + traceprod(dKxp_pseudo(:,:,i), Txp) ...
% $$$         + (c - V*KXxp*b)' * dKxp_pseudo(:,:,i) * b;
  end
end
