
% TESTBED_RECONSTRUCTIONS
%
% Plot the reconstruction.

% Copyright (c) 2010 Jaakko Luttinen

function hax = testbed_reconstructions(Q, varargin)

M = size(Q.W, 1);
N = size(Q.X, 2);

% Parse optional arguments
options = struct('stations',  1:M, ...
                 'times',     1:N, ...
                 'std_scale', 0);
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

% Reconstruction
Yh = bsxfun(@plus, Q.W(options.stations,:)*Q.X(:,options.times), ...
            Q.mu(options.stations));
if options.std_scale > 0
  varYh = zeros(size(Yh));
  for j=1:length(options.times)
    for i=1:length(options.stations)
      m = options.stations(i);
      n = options.times(j);
      varYh(i,j) = traceprod(Q.CovW(:,:,m),Q.CovX(:,:,n)) + Q.X(:,n)'* ...
          Q.CovW(:,:,m)*Q.X(:,n) + Q.W(m,:)*Q.CovX(:,:,n)*Q.W(m,:)' + Q.v_mu(m);
    end
  end
end

% Plot the reconstructions
time = Q.options.user_data.time(options.times);
if options.std_scale > 0
  hax = tserrorplot(time, Yh', 2*sqrt(varYh)');
else
  hax = tsplot(time, Yh);
end
N = length(hax);
set(hax, 'xlim', [min(time), max(time)]);
dateticks(hax, 'x', 'mmmyy', 'keeplimits');
set(hax(1:(N-1)), 'xticklabel', []);
set_subplot_positions(hax, N, 1, [0.05 0.05 0.05 0.1], [0.02 0.02]);

% Plot the data
addtsplot(time, Q.options.user_data.observations(options.stations,options.times), 'r');