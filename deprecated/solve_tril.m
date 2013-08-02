function X = solve_tril(L, Y)
% function X = solve_tril(L, Y)
%
% Solves X = L \ Y, or L*X = Y,  when L is lower triangular.

if issparse(L)
% $$$   if ~issparse(L)
% $$$     error('hei')
% $$$   end
  X = L \ Y;
else
  opts.LT = true;
  X = linsolve(L, full(Y), opts);
end
