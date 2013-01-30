
function aistats2012_run_comparison(p)

if nargin < 1 || p == 0
  p = 1:4;
end

% Noise std
s = 0.3;

%
% Theta
%

% Cholesky

if any(p==1)
  fid = fopen('/home/jluttine/matlab/publications/aistats2012/res1.txt', 'a');
% $$$   aistats2012_comparison(3, s, true, false, fid);
% $$$   aistats2012_comparison(3.5, s, true, false, fid);
% $$$   aistats2012_comparison(4, s, true, false, fid);
% $$$   aistats2012_comparison(4.5, s, true, false, fid);
  aistats2012_comparison(5, s, true, false, fid);
  fclose(fid);
end

% Kronecker

if any(p==2)
  fid = fopen('/home/jluttine/matlab/publications/aistats2012/res2.txt', 'a');
% $$$   aistats2012_comparison(3, s, true, true, fid);
% $$$   aistats2012_comparison(3.5, s, true, true, fid);
% $$$   aistats2012_comparison(4, s, true, true, fid);
% $$$   aistats2012_comparison(4.5, s, true, true, fid);
% $$$   aistats2012_comparison(5, s, true, true, fid);
% $$$   aistats2012_comparison(6, s, true, true, fid);
% $$$   aistats2012_comparison(7, s, true, true, fid);
  fclose(fid);
end

%
% F
%

% Kronecker

if any(p==3)
  fid = fopen('/home/jluttine/matlab/publications/aistats2012/res3.txt', 'a');
% $$$   aistats2012_comparison(3, s, false, true, fid);
% $$$   aistats2012_comparison(3.5, s, false, true, fid);
% $$$   aistats2012_comparison(4, s, false, true, fid);
% $$$   aistats2012_comparison(4.5, s, false, true, fid);
% $$$   aistats2012_comparison(5, s, false, true, fid);
% $$$   aistats2012_comparison(5.5, s, false, true, fid);
% $$$   aistats2012_comparison(6, s, false, true, fid);
% $$$   aistats2012_comparison(7, s, false, true, fid);
  aistats2012_comparison(8, s, false, true, fid);
  fclose(fid);
end

% Cholesky

if any(p==4)
  fid = fopen('/home/jluttine/matlab/publications/aistats2012/res4.txt', 'a');
% $$$   aistats2012_comparison(3, s, false, false, fid);
% $$$   aistats2012_comparison(3.5, s, false, false, fid);
% $$$   aistats2012_comparison(4, s, false, false, fid);
% $$$   aistats2012_comparison(4.5, s, false, false, fid);
% $$$   aistats2012_comparison(5, s, false, false, fid);
  aistats2012_comparison(5.5, s, false, false, fid);
  fclose(fid);
end
