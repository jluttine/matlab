% WSUM - Weighted sum of a matrix along one dimension.
%
% Y = WSUM(X,W,DIM)
%
% X   : matrix
% W   : vector of weights
% DIM : dimension to sum along (default: first non-singleton dimension)
%
% When X is a 2D-matrix, one can simply use W'*X or X*W for more
% efficient performace.

% Last modified 2010-06-09
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function y = wsum(X,w,dim)

% warning('This function seems to be inefficient..')

if ~isvector(w)
  error('Weights must be given as a vector.');
end

if nargin < 3
  % Default dimension: first non-singleton dimension
  dim = find(size(X)>1, 1);
  if isempty(dim)
    dim = 1;
  end
end

if size(X,dim) ~= length(w)
  error(['The length of the weight vector does not match with the size ' ...
         'of the matrix.']);
end

s = ones(1,max(dim,2));
s(dim) = length(w);
y = sum(bsxfun(@times, X, reshape(w,s)), dim);