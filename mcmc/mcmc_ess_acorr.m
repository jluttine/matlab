
% Effective sample size (ESS) based on autocorrelation
%
% Kass et al (1998), no, not this one

function Nh = mcmc_ess_acorr(x)

x = x(:);

N = length(x);
xh = mean(x);

g = zeros(size(x));
for k=0:(N-1)
%  g(k+1) = 1/(N-k) * (x(1:(N-k))-xh)'*(x((1+k):end));
  g(k+1) = 1/N * (x(1:(N-k))-xh)'*(x((1+k):end)-xh);
end


% $$$ % Cut-off
% $$$ rho = g/g(1);
% $$$ d = find(rho<0.01,1);
% $$$ tau = sum(rho(1:d));
% $$$ Nh = N / tau;
% $$$ figure
% $$$ plot(rho)

% $$$ figure
% $$$ plot(g/g(1))

%g_even = g(1:2:(end-1));
%g_odd = g(2:2:end);
%d = find(g_even+g_odd<=0,1)

d = find((g/g(1))<0.01, 1);
if isempty(d)
  d = length(g);
end

% $$$ g_sum = g(1:(end-1)) + g(2:end);
% $$$ d = find(g_sum<=0, 1);
% $$$ if isempty(d)
% $$$   d = length(g);
% $$$ end

%d
%v = 1/N * (g(1) + 2*sum(g(2:(2*d+1))));
v = 1/N * (g(1) + 2*sum(g(2:(d-1))));
Nh = g(1) / v;
%[~,d] = find(


% error('not implemented')