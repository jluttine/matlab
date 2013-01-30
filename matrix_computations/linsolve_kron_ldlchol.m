% LINSOLVE_KRON_LDLCHOL - Solves a matrix-vector equation when the matrix is
%                         a Kronecker product of two sparse symmetric
%                         positive definite matrices.
%
% Solves KRON(K1,K2)*VEC(Y) = VEC(X). Usage:
%
%   Y = LINSOLVE_KRON_LDLCHOL(LD1,LD2,X)
%
% where LD1 = LDLCHOL(K1) and LD2 = LDLCHOL(K2).
%
%   Y = LINSOLVE_KRON_LDLCHOL(...,FORM)
%
% FORM must be either 'vector' or 'matrix'. It specifies in which form the
% resulting Y is given. The default is 'matrix'.
%
% See also LDLCHOL, KRON.

% Last modified 2010-10-28
% Copyright (c) Jaakko Luttinen (jaakko.luttinen@tkk.fi)

function x = linsolve_kron_ldlchol(LD1, LD2, x, form);

K2 = @(x) ldlsolve(LD2,x);
K1 = @(x) ldlsolve(LD1,x);

if nargin < 4
  x = linsolve_kron(K1,K2,x);
else
  x = linsolve_kron(K1,K2,x,form);
end
