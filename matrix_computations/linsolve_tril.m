function X = linsolve_tril(L,Y)
% function X = LINSOLVE_TRIL(L, Y)
%
% Solves X = L \ Y, or L*X = Y,  when L is upper triangular.

% (c) 2010 Jaakko Luttinen

if issparse(L)
  X = L \ Y;
else
  opts.LT = true;
  X = linsolve(L, full(Y), opts);
end
