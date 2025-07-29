%% Figure 24. %% random realizations from covariance matrix
addpath ../data
addpath ../tools

load GOCE_NEATLANTIC
load GSHHS_i
load vik


minlon=-60;
maxlon= 30;
minlat= 55;
maxlat= 82.5;
%
[ny,nx,nr] = size(GOCE_NATL.m);
Nkl = numel(GOCE_NATL.m);
LL = chol(GOCE_NATL.Cm);
LL = LL';
sigma_rho=20^2;
%
figure,
for i=1:12
    
    subplot(4,3,i)
    Q_ksi = randn(Nkl,1);
    ksi = LL*Q_ksi;
    
    T=GOCE_NATL.m/std(GOCE_NATL.m(:))*sqrt(sigma_rho) + reshape(ksi,[ny nx nr]);
    
    plot_spher2(T,GOCE_NATL.lat,GOCE_NATL.lon,GOCE_NATL.r), colormap(flipud(vik))
    
    ishore = find([shorelines.Lon] > minlon &...
        [shorelines.Lon] < maxlon & ...
        [shorelines.Lat] > minlat & ...
        [shorelines.Lat] < maxlat);
    sh_lon = [shorelines.Lon];
    sh_lat = [shorelines.Lat];
    [sh_x,sh_y,sh_z]=sph2cart(sh_lon(ishore)*pi/180,sh_lat(ishore)*pi/180,sh_lat(ishore)*0+6375e3);
    plot3(sh_x,sh_y,sh_z,'.','Color',[.7 .7 .7])
    view(90,10)
    camlight('headlight'),
    drawnow
end

set(gcf,'Units','normalized','OuterPosition',[0 0 1 1])

print('-dpng','-r400','../fig/fig24')
print('-depsc','-r400','../fig/fig24')