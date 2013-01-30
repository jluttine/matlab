
function data = testbed_plot_data(data)

if nargin < 1
  data = testbed_loaddata();
  data = testbed_preprocess(data);
end

M = size(data.observations,1);
group = 10;
for m=1:group:M
  stations = m:min(m+group-1, M);
  tsplot(data.time, data.observations(stations, :), 'r');
end

if nargout == 0
  clear data;
end