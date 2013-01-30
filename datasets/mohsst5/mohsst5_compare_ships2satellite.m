function mohsst5_compare_ships2satellite

data = mohsst5_loaddata;
verdata = load('/home/alexilin/matlab/kaplan/verifdata', 'vd', 'lat', ...
               'lon', 'time');

i0 = find(data.time == verdata.time(1));
i1 = find(data.time == verdata.time(end));

time = data.time(i0:i1);
Y = data.observations(:,i0:i1);
Yver = verdata.vd;

figure
plot(Y(:), Yver(:), '.');
xlabel('Ship')
ylabel('Satellite')
set(gca, 'DataAspectRatio', [1 1 1]);
set(gca, 'XLimMode', 'manual');
set(gca, 'YLimMode', 'manual');

% Draw a line for comparison
hold on
x = [-100, 100];
y = [-100, 100];
plot(x(:), y(:), 'r');

ind = ~isnan(Y(:)) & ~isnan(Yver(:));
correlation = corr(Y(ind), Yver(ind))
Cov = cov(Y(ind)+10, Yver(ind))

eigenvalues = eig(Cov)
sqrt_eigenvalues = sqrt(eigenvalues)

inds = [200 300 299 750 1050 1270 1580];
hax = tsplot(time, Y(inds,:), 'r.-');
hold on
addtsplot(time, Yver(inds,:), 'k.-');
for i=1:length(hax)
  datetick(hax(i), 'x')
end