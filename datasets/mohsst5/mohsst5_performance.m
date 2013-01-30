
function mohsst5_performance(data, Qgp, Qpca)

Xr = Qgp.W * Qgp.X;
Vr = Qgp.W*Qgp.varX + Qgp.varW*Qgp.X + Qgp.varW*Qgp.varX;
lat = data.latitude;
lon = data.longitude;

file01 = '/home/jluttine/matlab/mohsst5/temp_perf01.mat';
save(file01, 'Xr', 'Vr', 'lat', 'lon');

if nargin < 3
  Qpca = load('/home/alexilin/matlab/kaplan/vbresults80');
end
Xr = bsxfun(@plus, Qpca.A * Qpca.S, Qpca.Mu);
Vr = zeros(size(Xr));

file02 = '/home/jluttine/matlab/mohsst5/temp_perf02.mat';
save(file02, 'Xr', 'Vr', 'lat', 'lon');

compare(file01, file02);

function compare(file1, file2)

%addpath(genpath('/home/alexilin/matlab'));

%clear all
load('/home/alexilin/matlab/kaplan/verifdata', 'vd', 'lat', 'lon', 'time');
time_vd = time;

load('/home/alexilin/matlab/kaplan/ssta', 'x', 'time');
load('/home/alexilin/matlab/kaplan/kaplanmask', 'M');
M1 = M;
M1(M1==0) = NaN;

recfile = {file1, file2};
% $$$ recfile = { 'vb80_rec', 'pcaoe_rec' };
% $$$ %recfile = { 'vb80_rec', 'vbd200_rec' };

warning off MATLAB:divideByZero
for i = 1:2
    load( recfile{i}, 'Xr', 'Vr', 'lat', 'lon' )
    [ctime,IA,IB] = intersect( time, time_vd );
    
    X2 = vd(:,IB);
    I{1} = ~isnan(x(:,IA)) & ~isnan(X2); % Part of train data
    I{2} = isnan(x(:,IA)) & ~isnan(X2);  % Test data
    I{3} = ~isnan(X2);                   % All data

    Ik{1} = ~isnan(x(:,IA)) & ~isnan(X2) & ...
            ~isnan(repmat(M1,1,length(IA))); % Part of train data
    Ik{2} = isnan(x(:,IA)) & ~isnan(X2) & ...
            ~isnan(repmat(M1,1,length(IA))); % Test data
    Ik{3} = ~isnan(X2) & ...
            ~isnan(repmat(M1,1,length(IA))); % All data
    
    [ sum(I{1}(:)) sum(I{2}(:)) sum(I{3}(:)) ...
      sum(Ik{1}(:)) sum(Ik{2}(:)) sum(Ik{3}(:)) ]
    
    fprintf( '%s\n', recfile{i} )
    
    W = sqrt(cos(lat/180*pi));
    for k = 1:length(I)
        X1 = Xr(:,IA);
        X1(~I{k}) = NaN; % Mark as NaN excluded parts of data

        err1{i,k} = (X1 - X2);
        err{i,k} = (X1 - X2) .* repmat(W,length(lon),size(X1,2));
        
        tmp = err{i,k}(:);
        N = sum( ~isnan(tmp) );
        tmp(isnan(tmp)) = 0;
        rms(i,k) = sqrt(sum(tmp.^2)/N);
        me(i,k) = sum(tmp)/N;
        stde(i,k) = sqrt(rms(i,k).^2 - me(i,k).^2);

        tmp = err1{i,k};
        N = sum( ~isnan(tmp), 2 );
        tmp(isnan(tmp)) = 0;
        rmsmap{i,k} = sqrt(sum(tmp.^2,2)./N);
        
        rmsts{i,k} = sqrt(sum(tmp.^2,1)./sum(~isnan(tmp),1));
        
        memap{i,k} = sum(tmp,2)./N;
        stdemap{i,k} = sqrt(rmsmap{i,k}.^2 - memap{i,k}.^2);
        
        tmp = stdemap{i,k};
        tmp(isnan(tmp)) = 0;
        stde_1(i,k) = sqrt( sum((tmp.^2).*N)/sum(N) );
        
        
        % Compute the same for Kaplan's areas
        err1_k{i,k} = err1{i,k} .* repmat(M1,1,size(X1,2));
        err_k{i,k} = err{i,k} .* repmat(M1,1,size(X1,2));
        
        tmp = err_k{i,k}(:);
        N = sum( ~isnan(tmp) );
        tmp(isnan(tmp)) = 0;
        rms_k(i,k) = sqrt(sum(tmp.^2)/N);
        me_k(i,k) = sum(tmp)/N;
        stde_k(i,k) = sqrt(rms_k(i,k).^2 - me_k(i,k).^2);

        tmp = err1_k{i,k};
        N = sum( ~isnan(tmp), 2 );
        tmp(isnan(tmp)) = 0;
        rmsmap_k{i,k} = sqrt(sum(tmp.^2,2)./N);
        
        memap_k{i,k} = sum(tmp,2)./N;
        stdemap_k{i,k} = sqrt(rmsmap_k{i,k}.^2 - memap_k{i,k}.^2);
        
        tmp = stdemap{i,k};
        tmp(isnan(tmp)) = 0;
        stde_1_k(i,k) = sqrt( sum((tmp.^2).*N)/sum(N) );

        switch k
            case 1, fprintf( 'Train data: ' )
            case 2, fprintf( 'Test data:  ' )
            case 3, fprintf( 'All data:   ' )
        end
        fprintf( '%i: rms=%.4f me=%.4f stde=%.4f stde1=%.4f ',...
                 i, rms(i,k), me(i,k), stde(i,k), stde_1(i,k) )
        fprintf( ' rms_k=%.4f me_k=%.4f stde_k=%.4f stde1_k=%.4f \n',...
                 rms_k(i,k), me_k(i,k), stde_k(i,k), stde_1_k(i,k) )
    end
    
    %wmapts( err{i}, lat, lon, ctime )
    %wmap( errmap{i}, lat, lon, 'clim', [ 0 max(errmap{i}(:))] )
    %colormap( jet )
end
warning on MATLAB:divideByZero

load obserr R
R = R(:,IA);

figure
plot( R(:), err{1,1}(:)-err{2,1}(:), '.' )

[Rsort,I1] = sort( R(:) );
d = err{1,1}(:)-err{2,1}(:);
d = d(I1);

hold on
plot( Rsort, filter( ones(1,100)/100, 1, d ), 'r' );

%return

if 0
    map = rmsmap{1,3};
    map(isnan(map)) = 0;
    wmap( map, lat, lon, 'clim', [ 0 2.5 ] )
    set( gcf, 'Name', 'RMS(VBPCA)' )
    colormap(1-gray)
    
    map = rmsmap{2,3};
    map(isnan(map)) = 0;
    wmap( map, lat, lon, 'clim', [ 0 2.5 ] )
    set( gcf, 'Name', 'RMS(PCAOE)' )
    colormap(1-gray)
end

d1 = rmsmap{1,3}-rmsmap{2,3};
d1(isnan(d1)) = 0;
wmap( d1, lat, lon, 'coolplot', 0 )
set( gcf, 'Name', 'RMS(VBPCA) - RMS(PCAOE)' )
figfontsize(7)

d1ts = rmsts{1,3}-rmsts{2,3};
plot_ts( d1ts, 'time', ctime, 'lstyle', 'k' )
set( gcf, 'Name', 'RMS(VBPCA) - RMS(PCAOE)' )
set( gcf, 'position', [ 360   478   560   125 ] )
set( gca, 'position', [ 0.13  0.15  0.775 0.80 ] )
figfontsize(12)

plot_ts( rmsts{2,3}, 'time', ctime, 'lstyle', 'k' )
set( gcf, 'Name', 'RMS(PCAOE)' )

dme = abs(memap{1,3})-abs(memap{2,3});
dme(isnan(dme)) = 0;
wmap( dme, lat, lon, 'coolplot', 0 )
set( gcf, 'Name', 'ME(VBPCA) - ME(PCAOE)' )

dstde = stdemap{1,3}-stdemap{2,3};
dstde(isnan(dstde)) = 0;
wmap( dstde, lat, lon, 'coolplot', 0 )
set( gcf, 'Name', 'STDE(VBPCA) - STDE(PCAOE)' )

mean(d1(~isnan(d1)))

return

d1m = d1;
d1m( M==0 ) = NaN;
wmap( d1m, lat, lon )
set( gcf, 'Name', 'VBPCA - PCAOE (mask)' )

err1m = err1;
err1m( M==0,: ) = NaN;
rms1m = sqrt(mean(err1m(~isnan(err1m)).^2));

err2m = err2;
err2m( M==0,: ) = NaN;
rms2m = sqrt(mean(err2m(~isnan(err2m)).^2));

rms1m
rms2m
