
function test_gradient_Kxx

mycheckgrad(@func, randn(10,1), 1e-8)


function [f, df] = func(x)

[K,tmp,dK_x] = gpcovDot(1,x',x');

f = sum(K(:));

n = length(x);
df = zeros(n,1);

Dpp = zeros(n,n);
for i=1:n
  %error('jou')
  Dpp(:) = 0;
  Dpp(i,:) = dK_x(1,:,i);
%  Dpp(:,i) = dK_x(1,:,i)';
  Dpp(:,i) = Dpp(:,i) + dK_x(1,:,i)';
%  Dpp(i,i) = 0; % DEBUGGING
  Dpp;
% $$$     Dpp(:,i) = dKpp_pseudo(j,:,i)';

    df(i) = sum(Dpp(:));
% $$$     traceprod(Dpp, Tpp, true) ...
% $$$         - 0.5 * b' * Dpp * b ...
% $$$         + dKxp_pseudo(j,:,i) * Txp(:,i) ...
% $$$         + dif' * dKxp_pseudo(j,:,i)' * b(i);
% $$$   end
end
