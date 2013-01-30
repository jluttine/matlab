% T_ML - ML estimation of the degrees of freedom of a multivariate Student-t
% distribution.
%
% Finds the maximum likelihood estimate for the degrees of freedom parameter
% of a (multivariate) Student-t probability density function
% STUDENT-T(Y|MU,COV,NU), where Y is the random variable, MU is the
% location, COV is the scale matrix and NU is the degrees of freedom.
%
% The function is called as:
%
% NU = T_ML(X2, NU_INIT, D)
%
% where
%
% X2      : (Y-MU)'*INV(COV)*(Y-MU)
% NU_INIT : the initial guess for NU
% D       : the dimensionality of the distribution
%
% NU : the maximum likelihood estimate

% Last modified 2010-06-08
% Copyright (c) Jaakko Luttinen

function nu = t_ml(x2, nu, D)

if nargin < 3
  D = 1;
end

%mycheckgrad(@cost, log(nu), 1e-4);

options = optimset('GradObj','on','Display','off','MaxIter',10);
nu = exp( fminunc(@cost, log(nu), options) );

% $$$ % DEBUG STUFF
% $$$ lognuh = -2:0.1:4;
% $$$ c = zeros(size(lognuh));
% $$$ dc = zeros(size(lognuh));
% $$$ for i=1:length(c)
% $$$   [c(i),dc(i)] = cost(lognuh(i));
% $$$ end
% $$$ figure(100)
% $$$ plot(lognuh, c)
% $$$ hold on
% $$$ plot(lognuh, dc, 'r')
% $$$ hold off
% $$$ 

  function [y, dy] = cost(lognu)
  nu = exp(lognu);
  [y,dy] = t_logpdf(x2, 0, nu, D);
  y = -sum(y);
  dy = -sum(dy)*nu;
  end
  
end