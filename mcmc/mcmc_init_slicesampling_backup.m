
% Axis aligned slice sampling

function func = mcmc_init_slicesampling(x_init, get_logpdf, varargin)

options = struct( ...
    'type',  'single-axis', ...  % single-axis / reflective-outside
    'y',     @(x) x,  ...        % transform the variable
    'logdy', @(x) 0,  ...        % log-jacobi of the transformation
    'dy',    @(x) 1,  ...
    'ddy',   @(x) 0,  ...
    'f',     @(y, varargin) []); % function for making pre-evaluations for logpdf

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);

% Initialization
x_current = x_init;
y_current = options.y(x_current);
fy_current = options.f(y_current);
logpdf_current = get_logpdf(y_current, fy_current);

x_proposal = x_current;

switch options.type
 
 case 'single-axis'
  func = @single_axis;
  % Initial width of the proposal
  w = ones(size(x_init));
  
 case 'reflective-outside'
  func = @reflective_outside;
  
 case 'reflective-inside'
  func = @reflective_inside;
  
 case default
  error('Unknown type for slice sampling');
  
end

  function [y, fy] = reflective_inside(varargin)
  
  end

  function [y, fy] = reflective_outside(varargin)
  disp('Check gradient')
  mycheckgrad(@helper, x_current, 1e-4)
  error('jou')
    function [f, df] = helper(x)
    y = options.y(x);
    [fy, dfy] = options.f(y);
    logpdf = get_logpdf(y, fy, dfy);
    [f, df] = logpdf(varargin{:});
    f = f + options.logdy(x);
    df = df.*options.dy(x) + options.ddy(x)./options.dy(x); % fix this!!
    end
  end

  function [y, fy] = single_axis(varargin)
  lp_current = logpdf_current(varargin{:}) + options.logdy(x_current);
  for i=1:numel(x_current)

    % Uniform vertical component
    lp_vertical = log(rand()) + lp_current;
    
    % DEBUGGING STUFF, DRAW LOGPDF
% $$$     if i == 2
% $$$       x1 = -10:0.1:5;
% $$$       x = repmat(x_current, [1, length(x1)]);
% $$$       x(i,:) = x1;
% $$$       lp = zeros(size(x1));
% $$$       for k=1:size(x,2)
% $$$         y = options.y(x(:,k));
% $$$         fy = options.f(y);
% $$$         logpdf = get_logpdf(y, fy);
% $$$         lp(k) = logpdf(varargin{:}) + options.logdy(x(:,k));
% $$$       end
% $$$       semilogx(exp(x(i,:)), lp);
% $$$       figure()
% $$$       imagesc(varargin{1})
% $$$       error('Drawing logpdf..')
% $$$     end
    
    % Find the initial horizontal interval
    j = 0;
    lp_min = inf;
    x_min = x_current;
    while lp_min >= lp_current;
      x_min(i) = x_current(i) - (j+1)*w(i);
      %x_min(i) = x_current(i) - 2^(j)*w(i);
      y_min = options.y(x_min);
      fy_min = options.f(y_min, i, y_current, fy_current);
      logpdf_min = get_logpdf(y_min, fy_min);
      lp_min = logpdf_min(varargin{:}) + options.logdy(x_min);
      j = j + 1;
      %disp('decrease min')
    end
    j = 0;
    lp_max = inf;
    x_max = x_current;
    while lp_max >= lp_current;
      x_max(i) = x_current(i) + (j+1)*w(i);
      %x_max(i) = x_current(i) + 2^(1*j)*w(i);
      y_max = options.y(x_max);
      fy_max = options.f(y_max, i, y_current, fy_current);
      logpdf_max = get_logpdf(y_max, fy_max);
      lp_max = logpdf_max(varargin{:}) + options.logdy(x_max);
      j = j + 1;
      %disp('increase max')
    end

    % Draw the horizontal component
    accepted = false;
    while ~accepted
      % Propose uniformly from the horizontal interval
      x_proposal(i) = unifrnd(x_min(i), x_max(i));
      y_proposal = options.y(x_proposal);
      % Get density value for the proposal
      fy_proposal = options.f(y_proposal, i, y_current, fy_current);
      logpdf_proposal = get_logpdf(y_proposal, fy_proposal);
      lp_proposal = logpdf_proposal(varargin{:}) + options.logdy(x_proposal);
      % Reject or accept the proposal
      if lp_proposal < lp_vertical
        % Reject the proposal and fix the horizontal interval
        %disp('reject')
        if x_proposal(i) < x_current(i)
          x_min(i) = x_proposal(i);
        else
          x_max(i) = x_proposal(i);
        end
      else
        % Accept the proposal
        logpdf_current = logpdf_proposal;
        lp_current = lp_proposal;
        x_current = x_proposal;
        y_current = y_proposal;
        fy_current = fy_proposal;
        accepted = true;
      end
    end
    % Store the width of the interval
    w(i) = (x_max(i) - x_min(i)) / 2;
  end
  %w
  y = y_current;
  fy = fy_current;
  end
end

