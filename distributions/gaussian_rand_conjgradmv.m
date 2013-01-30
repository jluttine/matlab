% GP_SAMPLE_PCG - Draw a sample from a conditional multivariate normal
%                 distribution using (preconditioned) conjugate gradient.
%
% The joint covariance matrix K is assumed to be of the form:
%
%   K = [  K1     K1
%          K1   (K1+K2) ]
%
% A sample is drawn from the conditional distribution
%
%   P(X1|X2) = N( X1 | K1*INV(K1+K2)*X2, K1 - K1*INV(K1+K2)*K1 )
%
% using (preconditioned) conjugate gradient method. This can be efficient
% for sparse covariance matrices.
%
% Usage:
%
% X1 = GP_SAMPLE_PCG(X2, K, K1, L1, L2)
% X1 = GP_SAMPLE_PCG(X2, K, K1, L1, L2, M)
%
% The matrices L1 and L2 are such that K1=L1*L1' and K2=L2*L2'. These can
% be found using, for instance, CHOL or LCHOL.
%
% One can also give function handles such that K(X), K1(X), L1(X) and L2(X)
% return the corresponding matrix-vector products.
%
% M is an optional preconditioner for the conjugate gradient method. To
% speed up the method, INV(M) should approximate INV(K1+K2). It can be
% extremely important to use a good preconditioner. One can also give a
% handle to a function which evaluates and returns M\X.
%
% Z1 = L1*RANDN(N1,1)
% Z2 = L2*RANDN(N2,1)
%
% One can also give a matrix S such that
%
%   K = [  K1       K1*S'
%         S*K1   S*(K1+K2)*S' ]
%
% and then a sample is drawn from the conditional Gaussian distribution
% with mean
%
%   E(X1|X2) = K1 * INV(S*(K1+K2)*S') * S*X2, 
%
% and covariance
%
%   COV(X1|X2) = K1 - K1 * S' * INV(S*(K1+K2)*S') * S * K1 )
%
% For instance, S can be an identity matrix from which some rows have been
% removed. This would correspond to having missing values in X2. In any
% case, S should be such that matrix-vector products are fast to evaluate,
% thus, preferably sparse. At least for now, S can not be a function
% handle. Note: the preconditioner might become worse when missing values..
%
% Usage:
%
% X1 = GP_SAMPLE_PCG(X2, K, K1, L1, L2, M, S)

% Last modified 2010-10-29
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function z1 = gaussian_rand_conjgradmv(K1, K2, z1, z2, I, varargin)
%function z1 = gaussian_rand_conjgradmv(x2, K, K1, z1, z2, I, varargin)
% To reduce memory usage: x1=z1
%
% Large variables in memory: x2, z1, z2

options = struct( ...
    'maxiter', 100,  ...
    'tol',     1e-6, ...
    'verbose', false);

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);

if nargin < 6 || isempty(I)
  I = true(size(x2));
end

if isnumeric(K1)
  K1 = @(x) K1*x;
end
if isnumeric(K2)
  K2 = @(x) K2*x;
end

% Assume zero means for now.
%mu1 = 0;
%mu2 = 0;

% Observations are subtracted from Z2!!
%z2 = z2 - x2;
%z2 = z2 - x2 + mu2
%z2 = z2 + z1 - x2 + mu2
z2 = conjgradmv(K2, z2, I, ...
                'tol', options.tol, ...
                'maxiter', options.maxiter, ...
                'verbose', options.verbose);
z1 = z1 - K1(z2);
%z1 = mu1 + z1 - K1(z2);

% x1 = z1;
