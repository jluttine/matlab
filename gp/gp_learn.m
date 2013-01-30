
% theta = gp_learn(y, covfunc, theta, varargin)

function [theta, loglikelihood] = gp_learn(y, covfunc, theta, varargin)

options = struct('maxiter', 100, ...
              'checkgrad', false);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

if options.checkgrad
  mycheckgrad(@cost, log(theta), 1e-3);
end

[logtheta, logpdfs] = minimize(log(theta), @cost, options.maxiter);
%options = optimset('GradObj', 'on');
%[logtheta, logpdfs] = fminunc(@cost, log(theta), options);

loglikelihood = -min(logpdfs);
theta = exp(logtheta);

  function [c, dc] = cost(logtheta)
  theta = exp(logtheta);
  %fprintf( '%f %f %f %f\n', theta(1), theta(2), theta(3), theta(4) )
  
  [c, dc] = gp_loglikelihood(y, covfunc, theta);
  c = -c;
  dc = -dc.*theta;
  
  if isinf(c)
      c = 1e0;
  end
%  c
  
  end

end
