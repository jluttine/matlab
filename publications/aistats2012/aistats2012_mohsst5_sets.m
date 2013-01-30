function [Ytrain, Ytest, Itrain, Itest, data, sea] = aistats2012_mohsst5_sets(testset, seed)

% Seed for the random number generator
if nargin < 2
  seed = 10;
end

rand('state', seed);
randn('state', seed);

% Load the data
data = mohsst5_loaddata();
sea = sum(~isnan(data.observations),2) > 0;
Y = data.observations(sea,:);

%
% FORM TRAIN AND TEST SETS
%

switch testset
 case 1 % randomly (uniformly) selected
  %disp('Test set: Uniform');
 
  I = rand(size(Y)) < 0.2;
  Itest = I & ~isnan(Y);
  Itrain = ~I & ~isnan(Y);
  Ytest = Y;
  Ytrain = Y;
  Ytest(~Itest) = nan;
  Ytrain(~Itrain) = nan;
  
  test_percent = sum(Itest(:)) / sum(~isnan(Y(:)));
  %fprintf('Size of the test set: %.1f percent\n', 100*test_percent);

 case 2 % pattern from earliest years used for latest years
  %disp('Test set: Pattern');
  
  Imv = isnan(Y);
  Imv(:,(20*12+1):end) = false;
  Itest = Imv(:,end:-1:1) & ~isnan(Y);
  Itrain = ~Imv(:,end:-1:1) & ~isnan(Y);
  
  Ytest = Y;
  Ytrain = Y;
  Ytest(~Itest) = nan;
  Ytrain(~Itrain) = nan;
  
  test_percent = sum(Itest(:)) / sum(~isnan(Y(:)));
  %fprintf('Size of the test set: %.1f percent\n', 100*test_percent);

 otherwise
  error('Unknown test set requested');
  
end
