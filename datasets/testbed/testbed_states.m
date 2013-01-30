
% TESTBED_TSPLOT
%
% Plot the latent signals.

% Copyright (c) 2010 Jaakko Luttinen

function hax = testbed_tsplot(Q, varargin)

% Parse optional arguments
options = struct();
[options, errmsg] = argparse(options, varargin{:});
error(errmsg);

% Plot the signals
time = Q.options.user_data.time(:);
hax = tserrorplot(time, Q.X', 2*sqrt(mdiag(Q.CovX))');
N = length(hax);
set(hax, 'xlim', [min(time), max(time)]);
dateticks(hax, 'x', 'mmmyy', 'keeplimits');
set(hax(1:(N-1)), 'xticklabel', []);
set_subplot_positions(hax, N, 1, [0.05 0.05 0.05 0.1], [0.02 0.02]);

if nargout < 1
  clear hax;
end
