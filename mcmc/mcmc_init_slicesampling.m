
% Axis aligned slice sampling

function [func, fx] = mcmc_init_slicesampling(x_init, get_logpdf, varargin)

options = struct( ...
    'fx',    @(x, varargin) x);

% Parse arguments
[options, errmsg] = argparse( options, varargin{:} );
error(errmsg);

% Initialization
x_current = x_init;
fx_current = options.fx(x_current);
logpdf_current = get_logpdf(fx_current);

func = @single_axis;
fx = fx_current;

x_proposal = x_current;

w = 0.5*ones(size(x_current));
min0 = x_current - w;
max0 = x_current + w;

  function [x, fx] = single_axis(varargin)
  lp_current = logpdf_current(varargin{:});
  for i=1:numel(x_current)

    % Uniform vertical component
    lp_vertical = log(rand()) + lp_current;
    
% $$$     % DEBUGGING STUFF, DRAW LOGPDF
% $$$     if i == 1
% $$$       x1 = -2:0.1:2;
% $$$       x = repmat(x_current, [1, length(x1)]);
% $$$       x(i,:) = x1;
% $$$       lp = zeros(size(x1));
% $$$       for k=1:size(x,2)
% $$$         [fx, logpdf] = get_logpdf(x(:,k));
% $$$         lp(k) = logpdf(varargin{:});
% $$$       end
% $$$       plot(x(i,:), lp);
% $$$       error('Drawing logpdf..')
% $$$     end
    
    if isinf(lp_current) || isnan(lp_current)
      lp_current
      error('What?!?! Has accepted Nan/Inf logpdf???')
    end
    
    % Find the initial horizontal interval
    j = 0;
    lp_min = inf;
    x_min = x_current;
    while lp_min >= lp_vertical && j < 10
      x_min(i) = min0(i) - j*w(i);
      fx_min = options.fx(x_min, i, fx_current);
      logpdf_min = get_logpdf(fx_min);
      lp_min = logpdf_min(varargin{:});
      j = j + 1;
      %disp('decrease min')
    end
    j = 0;
    lp_max = inf;
    x_max = x_current;
    while lp_max >= lp_vertical && j < 10
      x_max(i) = max0(i) + j*w(i);
      fx_max = options.fx(x_max, i, fx_current);
      logpdf_max = get_logpdf(fx_max);
      lp_max = logpdf_max(varargin{:});
      j = j + 1;
      %disp('increase max')
    end

    % Draw the horizontal component
    accepted = false;
    while ~accepted
      %[i, x_min(i), x_max(i)]
      % Propose uniformly from the horizontal interval
      x_proposal(i) = unifrnd(x_min(i), x_max(i));
      % Get density value for the proposal
      fx_proposal = options.fx(x_proposal, i, fx_current);
      logpdf_proposal = get_logpdf(fx_proposal);
      lp_proposal = logpdf_proposal(varargin{:});
      %[lp_proposal, lp_vertical, lp_current]
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
        if isinf(lp_proposal)
          x_proposal
          lp_proposal
          error('Debugging WHAAAT?!?! Got INF logpdf..?!');
        end
        % Accept the proposal
        x_current = x_proposal;
        fx_current = fx_proposal;
        logpdf_current = logpdf_proposal;
        lp_current = lp_proposal;
        accepted = true;
      end
    end
    % Store the width of the interval
    min0(i) = x_min(i);
    max0(i) = x_max(i);
    w(i) = 0.5*(x_max(i) - x_min(i));
  end
  x = x_current;
  fx = fx_current;
  end

end
