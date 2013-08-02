function H = mvnentropy(Cov, ischol)
% H = entropy_gaussian(Cov)
%
% Cov is covariance matrix, or optionally, an upper triangular matrix U such
% that U'*U = Cov (Cholesky decomposition) which is marked ischol=true. The
% method calculates the Cholesky decomposition if it is not given.

% $$$ warning(['Function mvnentropy deprecated, start using gaussian_entropy (it uses' ...
% $$$          ' different parameters)'])

if nargin < 2 || ~ischol
  U = safechol(Cov, 1e-8);
else
  U = Cov;
end

d = size(Cov,1);

H = d/2*log(2*pi) + 0.5*logdetchol(U) + d/2;
