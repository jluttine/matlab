function Q = gpfa(D, Y, W_module, X_module, noise_module, varargin)

% OK, a dummy wrapper at the moment, but I could make this method a bit
% easier to use. For instance, this function would form the modules given
% some covariance functions or something...

Q = vbfa(D, Y, W_module, X_module, noise_module, varargin{:});
