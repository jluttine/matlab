function mohsst5_reconstruction_movie(Yh, errYh)

warning('Under construction')

data = mohsst5_loaddata();
Y = data.observations;

N = size(Yh,1);

c1 = max(abs(Y(:)));
c2 = max(abs(Yh(:)));
cl = max(c1,c2);
cl = prctile(abs([Yh(:);Y(:)]), 99.9)

cl_err = max(abs(errYh(:)));

N = 4;
for n=1:N
  figure
  hax = subplot(3,1,1);
  mohsst5_mapplot(Y(:,n), 'colormap', 'centered');
  set(hax, 'clim', [-cl cl]);
  hax = subplot(3,1,2);
  mohsst5_mapplot(Yh(:,n), 'colormap', 'centered');
  set(hax, 'clim', [-cl cl]);
  hax = subplot(3,1,3);
  mohsst5_mapplot(errYh(:,n), 'colormap', 'centered');
  set(hax, 'clim', [-cl_err cl_err]);
end