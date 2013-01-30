
function gpmapcolor(coord, f, Covf, lon, lat, logtheta, covfunc)

if iscell(coord)
  D = numel(coord);
end

%[lon, lat] = lonlat2standard(lon, lat);
[LAT,LON] = meshgrid(lat,lon);
Xh = [LON(:)'; LAT(:)'];

for d=1:D 
  subplot(D,2,2*d-1);
  %mapproj('global-ellipse');
  [fh, varfh] = gppred(coord{d}, f{d}, Covf{d}, Xh, logtheta{d}, covfunc{d});
  m_pcolor(LON,LAT,reshape(fh,size(LON)));
  shading('flat');
  mapcoast
  %colormap('hsv');

  subplot(D,2,2*d);
  %mapproj('global-ellipse');
  [fh, varfh] = gppred(coord{d}, f{d}, Covf{d}, Xh, logtheta{d}, covfunc{d});
  m_pcolor(LON,LAT,reshape(varfh,size(LON)));
  shading('flat');
  mapcoast
  %colormap('gray');
end

