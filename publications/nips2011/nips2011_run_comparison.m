
function aistats2012_run_comparison

%
% Theta
%

% Cholesky

nips2011_comparison(3, 1, true, false);
nips2011_comparison(3.5, 1, true, false);
nips2011_comparison(4, 1, true, false);
nips2011_comparison(4.5, 1, true, false);
nips2011_comparison(5, 1, true, false);

% Kronecker

% $$$ nips2011_comparison(3, 1, true, true);
% $$$ nips2011_comparison(3.5, 1, true, true);
% $$$ nips2011_comparison(4, 1, true, true);
% $$$ nips2011_comparison(4.5, 1, true, true);
% $$$ nips2011_comparison(5, 1, true, true);
% $$$ nips2011_comparison(6, 1, true, true);
nips2011_comparison(7, 1, true, true);

%
% F
%

% Kronecker

% $$$ nips2011_comparison(3, 1, false, true);
% $$$ nips2011_comparison(3.5, 1, false, true);
% $$$ nips2011_comparison(4, 1, false, true);
% $$$ nips2011_comparison(4.5, 1, false, true);
% $$$ nips2011_comparison(5, 1, false, true);
% $$$ nips2011_comparison(5.5, 1, false, true);
% $$$ nips2011_comparison(6, 1, false, true);
% $$$ nips2011_comparison(7, 1, false, true);
% $$$ nips2011_comparison(8, 1, false, true);

% Cholesky

nips2011_comparison(3, 1, false, false);
nips2011_comparison(3.5, 1, false, false);
nips2011_comparison(4, 1, false, false);
nips2011_comparison(4.5, 1, false, false);
nips2011_comparison(5, 1, false, false);
nips2011_comparison(5.5, 1, false, false);

