
function test_gp_geographic

lat = -87.5:5:87.5;
lon = -177.5:5:177.5;


%[LAT,LON] = meshgrid(lat,lon);
[LON,LAT] = meshgrid(lon,lat);
X = geographic_to_euclidean([LON(:)';LAT(:)']);


switch 1
 case 1
  % Use block-Toeplitz structure
  [lat,lon0] = meshgrid(lat,lon(1));
  X0 = geographic_to_euclidean([lon0(:)';lat(:)']);
  D = sqrt(sq_dist(X0,X));
  covfunc = gp_cov_toeplitz_block(gp_cov_pp(D,3));
 case 2
  % Brute force
  D = sqrt(sq_dist(X,X));
  covfunc = gp_cov_pp(D,3);
end

% covfunc1 = gp_cov_scale(covfunc1);
%covfunc1 = gp_cov_toeplitz(covfunc1);

theta = 4000;

t = cputime();
K = covfunc(theta);
time_covfunc = cputime() - t

figure
%spy(K)

fullness_K = nnz(K) / numel(K)
%return

t = cputime();
[L,p,q] = lchol(K);
time_lchol = cputime() - t

t = cputime();
y(q) = L*randn(size(K,1),1);
time_rand = cputime() - t

N = numel(y)

%figure
map_projection()
%map_pcolor(lon, lat, reshape(K(415,:), [length(lon), length(lat)]))
map_pcolor(LON, LAT, reshape(y, size(LON)))
%map_pcolor(lon, lat, reshape(y, [length(lon), length(lat)]))
map_colormap()
map_grid()
map_coast()
