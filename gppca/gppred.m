%
% [fh, varfh] = gppred(X, f, Covf, Xh, logtheta, covfunc, varargin)
%
function [fh, varfh] = gppred(X, f, Covf, Xh, logtheta, covfunc, ...
                              varargin)

warning('Deprecated. Use gp_predict or gp_predict_pseudo.');

% If you need more than just variances for PCA, you could give some sparse
% matrix to indicate which covariances should be evaluated. But again, this
% covariances should be zero for the latent function values fh a
% priori.. Things would be much easier if I just factored the principal
% components too.. :)

if ~iscell(covfunc)
  covfunc = {covfunc};
end

opts = struct('cholkx', [], 'khx', [], 'kh', []);

[ opts, errmsg, wrnmsg ] = argschk( opts, varargin{:} );
if ~isempty(errmsg), error( errmsg ), end
if ~isempty(wrnmsg), warning( wrnmsg ), end

nh = cols(Xh);

I = diag(sparse(ones(cols(X),1)));

if isempty(opts.cholkx)
  Kpp = feval(covfunc{:}, logtheta, X, X);

  % DEBUG: REGULARIZE
  %Kpp = regularize(Kpp);

  Kpp = regularize(Kpp);
  Lp = chol(Kpp, 'lower');
%  [Lp,Kpp] = safechol(Kpp, 1e-10, 'lower');
  % DEBUG:
  %Kpp = Lp*Lp';
else
  Lp = opts.cholkx;
end


if isempty(opts.khx)
  Kxp = feval(covfunc{:}, logtheta, Xh, X);
else
  Kxp = opts.khx;
end

if isempty(opts.kh)
  Kx = feval(covfunc{:}, logtheta, Xh, []);
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

% $$$ if isequal(covfunc{1}, @gpcovPP)
% $$$ %if strcmpi(covfunc{1}, @gpcovPP)
% $$$   if ~issparse(Lp)
% $$$     error('Lp should be sparse');
% $$$   else
% $$$     density_in_gppred = nnz(Lp) / prod(size(Lp))
% $$$   end
% $$$ end

