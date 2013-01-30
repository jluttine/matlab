function X = regularize(X, reg)

%error('Do not use this function');

% $$$ normX = normest(X, 1e-1);

n = length(X);
I = diag(sparse(ones(n,1)));

if nargin < 2
  reg = 1;
end
normX = reg;
coeff = 1e-6;
X = X + coeff*normX*I;

% $$$ while ~success
% $$$   try
% $$$     U = chol(X + 1e-8*normX*I, varargin{:});
% $$$     success = true;
% $$$   catch
% $$$     times = times + 1;
% $$$     warning(sprintf(['Cholesky decomposition failed %d-th time, increasing ' ...
% $$$                      'regularization'], times))
% $$$     normX = normX*100;
% $$$   end
% $$$ end
