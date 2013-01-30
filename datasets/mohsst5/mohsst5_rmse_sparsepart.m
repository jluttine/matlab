function weighted_rmse = mohsst5_rmse_sparsepart(Yh, data)

Itest = load('ind20test.mat');
Itest = Itest.Itest;
if nargin < 2 || isempty(data)
  data = mohsst5_loaddata;
end

% Test set
Ytest = data.observations;
Ytest(~Itest) = nan;
missing = sum(isnan(Ytest(:))) / numel(Ytest);

% Take only data from EARLY parts
rows = size(data.observations,1);
rminds = data.time > datenum([1900 1 1 0 0 0]); % remove by time period
rminds = rowsum(~isnan(data.observations)) > 0.3 * rows; % remove by sparsity
rminds = data.time > datenum([9999 1 1 0 0 0]); % full data
Ytest(:,rminds) = nan;
missing = sum(isnan(Ytest(:))) / numel(Ytest);
testsize = sum(~isnan(Ytest(:)));

% $$$ % Take only data from sparse parts
% $$$ rows = size(data.observations,1);
% $$$ rminds = rowsum(~isnan(data.observations)) > 0.7 * rows;
% $$$ Ytest(:,rminds) = nan;
% $$$ missing = sum(isnan(Ytest(:))) / numel(Ytest)
% $$$ testsize = sum(~isnan(Ytest(:)))

weights = sqrt(cosd(data.coordinates(2,:)));
weights = weights(:);
Weights = bsxfun(@times, ~isnan(Ytest), weights(:)); % mask with missing values
Yhw = bsxfun(@times, Yh, weights(:));
Ytestw = bsxfun(@times, Ytest, weights(:));
err = Yhw(:) - Ytestw(:);
err(isnan(err)) = []; % remove missing values
weighted_rmse = sqrt(err'*err / sum(Weights(:).^2));
