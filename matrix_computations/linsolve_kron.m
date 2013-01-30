% LINSOLVE_KRON - Solves a matrix-vector equation when the matrix is a
%                 Kronecker product of two matrices.
%
% Solves KRON(K1,K2)*VEC(Y) = VEC(X). Usage:
%
%   Y = LINSOLVE_KRON(K1,K2,X)
%
% K1 and/or K2 can be handles to functions which return the solutions K1\Z
% and K2\Z for given matrix inputs Z. By default, MLDIVIDE is used.
%
%   Y = LINSOLVE_KRON(...,FORM)
%
% FORM must be either 'vector' or 'matrix'. It specifies in which form the
% resulting Y is given. The default is 'matrix'.
%
% Y must be in matrix form.
%
% See also KRON.

% Last modified 2010-10-28
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function x = linsolve_kron(K1, K2, x, form);

% x = reshape(x,size(K2,1),size(K1,1));

if isnumeric(K2)
  x = K2\x;
else
  x = K2(x);
end
if isnumeric(K1)
  x = (K1\(x'))';
else
  x = K1(x')';
end

if nargin >= 4 && ~isempty(form)
  if strcmpi(form,'vector')
    x = x(:);
  elseif ~strcmpi(form,'matrix')
    error(['Requested the result in unknown form. Must be either ''vector'' ' ...
           'or ''matrix''']);
  end
end
