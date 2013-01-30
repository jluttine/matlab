%
% [fh, varfh] = gppred(X, f, Covf, Xh, kX, logtheta, ...)
%
% p(y|X) = N(y; f(X), inv(V))
function [f, varf] = spgpr(logtheta, covfunc, X, y, V, varargin)
%function [fh, varfh] = gppred(X, f, Covf, Xh, kX, logtheta, varargin)

% If you need more than just variances for PCA, you could give some sparse
% matrix to indicate which covariances should be evaluated. But again, this
% covariances should be zero for the latent function values fh a
% priori.. Things would be much easier if I just factored the principal
% components too.. :)

opts = struct('cholkx', [], 'khx', [], 'kh', []);
%opts = struct('cholkx', [], 'khx', [], 'kh', []);

[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end

nh = cols(Xh);

if isempty(opts.cholkx)
  Kpp = feval(kX, X, X, logtheta);
  Kpp = Kpp + 1e-6*eye(size(Kpp));
  Lp = chol(Kpp, 'lower');
else
  Lp = opts.cholkx;
end

if isempty(opts.khx)
  Kxp = feval(kX, Xh, X, logtheta);
else
  Kxp = opts.khx;
end

if isempty(opts.kh)
  Kx = feval(kX, Xh, [], logtheta);
else
  Kx = opts.kh;
end

fh = Kxp * solve_triu(Lp', solve_tril(Lp, f));

if nargout >= 2
  varfh = zeros(nh,1);
  for i=1:nh
    r = solve_tril(Lp, Kxp(i,:)');
    s = solve_triu(Lp', r);
    varfh(i) = Kx(i) - r'*r + s'*Covf*s;
  end
end

