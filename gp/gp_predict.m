
% [f_mean, f_var] = gp_pred(y, K_y, K_yf, k_f)

function [f_mean, f_var] = gp_predict(y, K_y, K_yf, k_f)

if issparse(K_y)
% $$$   LD = ldlchol(K_y);
% $$$   linsolve = @(x) linsolve_ldlchol(LD, x);
  % Don't use q, because linsolve(K_yf) would take a lot of time for
  % permutations..
% $$$   [L,p,q] = lchol(K_y);
% $$$   if p~=0
% $$$     p
% $$$     error('Matrix K_y not positive definite.');
% $$$   end
% $$$   linsolve = @(x) linsolve_lchol(L, x, q);
  [LD,p,q] = ldlchol(K_y);
  if p~=0
    p
    error('Matrix K_y not positive definite.');
  end
  linsolve = @(x) linsolve_ldlchol(LD, x, q);
else
  L = chol(K_y, 'lower');
  linsolve = @(x) linsolve_lchol(L, x);
end

f_mean = K_yf' * linsolve(y);
if nargout >= 2
  if issparse(K_yf)
% $$$     f_var = zeros(size(k_f));
% $$$     invK_y = spinv(K_y);
% $$$     for n=1:size(K_yf,2)
% $$$       f_var(n) = k_f(n) - K_yf(:,n)'*(invK_y*K_yf(:,n));
% $$$     end
    for n=1:size(K_yf,2)
      f_var(n) = k_f(n) - K_yf(:,n)'*linsolve(K_yf(:,n));
    end
  else
    f_var = k_f - dot(K_yf, linsolve(K_yf), 1)';
  end
end
  
% $$$   K_y = K_f + V;
% $$$   
% $$$   if issparse(K_y)
% $$$     LD = ldlchol(K_y);
% $$$     linsolve = @(x) linsolve_ldlchol(LD, x);
% $$$   else
% $$$     L = chol(K_y, 'lower');
% $$$     linsolve = @(x) linsolve_lchol(L, x);
% $$$   end
% $$$   
% $$$   f_mean = K_f * linsolve(y);
% $$$   if nargout >= 2
% $$$     f_var = diag(K_f) - dot(K_f, linsolve(K_f), 1)';
% $$$   end
  
