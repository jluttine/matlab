function [F,FF] = icml2011_mohsst5_parse_results(filename, N, burnin)

F_mean = 0;
FF = 0;

for n=(burnin+1):N
  fprintf('%d/%d\n', n, N);
  name = sprintf('%s_F%d',filename,n);
  load(name);
  F_mean = (F + (n-burnin-1)*F_mean) / (n-burnin);
  FF = (F.*F + (n-burnin-1)*FF) / (n-burnin);
end

F = F_mean;

