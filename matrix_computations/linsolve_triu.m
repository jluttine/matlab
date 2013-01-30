function X = linsolve_triu(R, Y, transpose)
% function X = LINSOLVE_TRIU(R, Y)
%
% Solves X = R \ Y, or R*X = Y,  when R is upper triangular.

% (c) 2010 Jaakko Luttinen

if issparse(R)
  X = R \ Y;
else
  opts.UT = true;
  if nargin >= 3 && transpose
    opts.UT = false;
    opts.LT = true;
    opts.TRANSA = true;
  end
  X = linsolve(R, full(Y), opts);
end
