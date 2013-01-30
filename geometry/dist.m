% D = DIST(X1, X2)
%
% Returns a matrix D of pair-wise distances between X1 and X2.
%
% X1 is K x M matrix
% X2 is K x N matrix
%
% K is the dimensionality, M and N are the number of samples.  The resulting
% matrix D is M x N matrix.
%
% Optional arguments:
%
% 'distfunc' : uses the given distance measure between the vectors
%              (default is the second norm).

% Copyright (c) 2010 Jaakko Luttinen

function D = dist(X1, X2, varargin)

% Default parameters
options = struct( ...
    'distfunc', @(z1,z2) norm(z1-z2)); % distance measure

% Check arguments
options = argget( options, varargin{:} );

% Evaluate pair-wise distances
n1 = size(X1,2);
n2 = size(X2,2);
D = zeros(n1,n2);
for i=1:n1
  for j=1:n2
    D(i,j) = options.distfunc(X1(:,i),X2(:,j));
  end
end
