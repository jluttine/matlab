
% Reflective multivariate slice sampling

function [func, fx] = mcmc_init_reflective(x_init, get_logpdf, get_dlogpdf, ...
                                           epsilon, L, varargin)

options = struct( ...
    'fx',    @func_x_default, ...
    'type',  'inside'); 
% $$$     'y',     @(x) x,  ...        % transform the variable
% $$$     'logdy', @(x) 0,  ...        % log-jacobi of the transformation
% $$$     'dy',    @(x) 1,  ...
% $$$     'ddy',   @(x) 0,  ...
% $$$     'f',     @(y, varargin) []); % function for making pre-evaluations for logpdf

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);

% Initialization
x_current = x_init;
[fx_current, dfx_current] = options.fx(x_current);
logpdf_current = get_logpdf(fx_current);

fx = fx_current;

x_proposal = x_current;

switch options.type
 
 case 'outside'
  func = @reflective_outside;
  
 case 'inside'
  func = @reflective_inside;
  
 case default
  error('Unknown type for slice sampling');
  
end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  function [x, fx] = reflective_inside(varargin)
  x_trajectory = x_current;
  fx_trajectory = fx_current;
  dfx_trajectory = dfx_current;
  logpdf_trajectory = logpdf_current;
  lp_trajectory = logpdf_current(varargin{:});

  % Momentum stuff
  step = exprnd(epsilon); % step length
  p = randn(size(x_current));
  
  %traj = [];

  % Simulate N steps
  for n=1:L
    
    if mod(n-1,ceil(L/10)) == 0
      % Update the slice
      lp_slice = lp_trajectory + log(rand());
    end
  
    % Move a step forward
    x_proposed = x_trajectory + step*p;
    % Check whether the new point is inside the slice
    [fx_proposed, dfx_proposed] = options.fx(x_proposed);
    logpdf_proposed = get_logpdf(fx_proposed);
    lp_proposed = logpdf_proposed(varargin{:});
    if lp_proposed <= lp_slice
      % Reflect based on the gradient at the inside point
      
      % Get the gradient
      dlogpdf = get_dlogpdf(dfx_trajectory);
      g = dlogpdf(varargin{:});

      % New momentum
      g2 = g'*g;
      if g2 > 0
        p = p - (2*(p'*g)/g2)*g;
      else
        % Random or opposite direction?
        p = -p;%orth(randn(numel(p))) * p;
      end

      % Check the reflected point is inside the slice
      x_reflect = x_trajectory + step*p;
      [fx_reflect, dfx_reflect] = options.fx(x_reflect);
      logpdf_reflect = get_logpdf(fx_reflect);
      lp_reflect = logpdf_reflect(varargin{:});
      if lp_reflect <= lp_slice
        disp('Reject trajectory')
        % TODO: Should I keep the state or ignore
        x = x_current;
        fx = fx_current;
        %x = nan;
        %fx = nan;
        return
      end
      
      %
      % TODO: HOW TO REQUEST THE GRADIENT?!?!?!!!!
      %

      % Verify opposite direction
      x_verify = x_trajectory - step*p;
      fx_verify = options.fx(x_verify);
      logpdf_verify = get_logpdf(fx_verify);
      lp_verify = logpdf_verify(varargin{:});
      if lp_verify > lp_slice
        disp('Reject trajectory')
        % TODO: Should I keep the state or ignore
        x = x_current;
        fx = fx_current;
        %x = nan;
        %fx = nan;
        return
      end
      
      % Accept the reflected step
      %disp('Reflect')
      x_trajectory = x_reflect;
      fx_trajectory = fx_reflect;
      dfx_trajectory = dfx_reflect;
      lp_trajectory = lp_reflect;
      logpdf_trajectory = logpdf_reflect;

    else
      
      % Accept the forward step
      x_trajectory = x_proposed;
      fx_trajectory = fx_proposed;
      dfx_trajectory = dfx_proposed;
      lp_trajectory = lp_proposed;
      logpdf_trajectory = logpdf_proposed;
      
    end
%    traj = [traj, x_trajectory];
  end

% $$$   figure
% $$$   plot(traj(1,:),traj(2,:), 'kx-');
% $$$   % Plot contours
% $$$   v = axis();
% $$$   vx = linspace(v(1), v(2), 20);
% $$$   vy = linspace(v(3), v(4), 10);
% $$$   [VY,VX] = meshgrid(vy,vx);
% $$$   V = [VX(:), VY(:)]';
% $$$   F = zeros(size(VX));
% $$$   for n=1:numel(F)
% $$$     logpdf = get_logpdf(V(:,n));
% $$$     F(n) = logpdf();
% $$$   end
% $$$   hold on
% $$$   contour(VX,VY,reshape(F,size(VX)));
% $$$   error('jou')

  disp('Accept trajectory')
  x_current = x_trajectory;
  fx_current = fx_trajectory;
  dfx_current = dfx_trajectory;
  logpdf_current = logpdf_trajectory;
  
  x = x_current;
  fx = fx_current;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function [y, fy] = reflective_outside(varargin)
  disp('Check gradient')
  mycheckgrad(@helper, x_current, 1e-6)
  error('Stopped because wanted to check gradients')
    function [f, df] = helper(x)
    [fx, dfx] = options.fx(x);
    logpdf = get_logpdf(fx);
    dlogpdf = get_dlogpdf(dfx);
    f = logpdf(varargin{:});
    df = dlogpdf(varargin{:});
    end
  end

end

function [fx, dfx] = func_x_default(x)
fx = x;
dfx = x;
end
