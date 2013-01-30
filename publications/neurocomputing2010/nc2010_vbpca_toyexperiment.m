
% Test whether it is necessary to rotate VB PCA
function nc2010_vbpca_toyexperiment(n,m,d,ncomps,rotate,maxiters, ...
                                    startupdatehyper,flatW,datatype)

if nargin < 3
  n = 200;
  m = 10;
  d = 4;
end

if nargin < 7
  startupdatehyper = 1;
end
if nargin < 6
  maxiters = 10000;
end
if nargin < 5
  rotate = true;
end
if nargin < 4
  ncomps = d;
end

randn('state', 1);
rand('state', 1);

% Generate multivariate normal data
mu = randn(m,1);
R = orth(randn(m,m)); % arbitrary rotation
switch datatype
 case 1
  disp('Using weak subspace')
  eigs = (1 + [d:-1:1 zeros(1,m-d)]) .^ 2
  datastring = 'weak';
 case 2
  disp('Using strong subspace')
  eigs = ([5*ones(1,d) 1*ones(1,m-d)]) .^ 2
  datastring = 'strong';
 case 3
  disp('Using no subspace')
  eigs = (m:-1:1) .^ 2
  datastring = 'no';
end  
Cov = R*diag(eigs)*R';
Yn = mvnrnd(mu', Cov, n)';

% Generate some outliers (NO!)
Yno = Yn;

% Generate some missing values
missing_values = 0.2
Ynom = Yno;
pmv = (rand(m,n) < missing_values);
Ytest = Ynom;
Ynom(pmv) = NaN;
Ytest(~pmv) = NaN;

if flatW
  disp('Using non-hierarchical flat prior for W');
  stringW = 'flatW';
else
  disp('Using hierarchical prior for W');
  stringW = 'hierW';
end


filename = sprintf(['/home/jluttine/matlab/neurocomputing2010/' ...
                    'nc2010_vbpca_toyexperiment_%s_subspace=%s_n=%d_m=%d_d=%d_ncomps=' ...
                    '%d_rotate=%d'], stringW,datastring,n,m,d,ncomps, rotate);


results = vbpcamv(Ynom, ncomps, 'testset', Ytest, 'maxiters', maxiters, ...
                  'rotate', rotate, 'startupdatehyper', startupdatehyper, ...
                  'startrotate', 1, 'autosavetime', 1, 'autosavefile', ...
                  filename, 'fixw', flatW);
  
