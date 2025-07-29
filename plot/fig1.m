%clear all 
%close all
%
addpath ../data
load GoceGrids% GOCE gradients from grids
load topography %topo
load GSHHS_i %coastlines
%
minlat = 55; maxlat = 82; nlat=maxlat-minlat;
minlon = -55; maxlon = 20; nlon=maxlon-minlon;
study_area = [linspace(minlon,maxlon,100)',maxlat*ones(100,1);...
    linspace(maxlon,minlon,100)',minlat*ones(100,1);...
    minlon,maxlat];
%
NoL = 1801;
NoP = 900; 
LAT = reshape(T.data(:,2),NoL,NoP);
LON = reshape(T.data(:,1),NoL,NoP);
% gradients
XXA = reshape(T.data(:,3),NoL,NoP);
XYA = reshape(T.data(:,4),NoL,NoP);
XZA = reshape(T.data(:,5),NoL,NoP);
ZZA = reshape(T.data(:,6),NoL,NoP);
YYA = reshape(T.data(:,8),NoL,NoP);
YZA = reshape(T.data(:,7),NoL,NoP);
%
ny0=1081; nx0=2161;
lontop = linspace(-180,180,nx0); lattop = linspace(  90,-90,ny0);
[Phi,Lam] = meshgrid(lattop,lontop);
topo_int = interp2(Phi,Lam,topoi',LAT, LON);
%
figure(1),clf
subplot(2, 2, 1);
axesm('ortho','origin',[64 -22]), axis off
surfm(LAT, LON, XXA, topo_int*1e-5-1); 
camlight, material dull, lighting Gouraud, shading interp, 
title('(a)            XX             '); 
c = colorbar('Position',[0.475641256296906 0.571251548946716 0.0117370892018775 0.151177199504337]); 
ylabel(c, 'E'); 
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp, 
%colormap(colbrew(1))
plotm(study_area(:,2),study_area(:,1),'b','LineWidth',1.5),
caxis([-1 1])
subplot(2, 2, 2);
axesm('ortho','origin',[64 -22]), axis off
surfm(LAT, LON, YYA, topo_int*1e-5-1); 
camlight, material dull, lighting Gouraud, shading interp
title('(b)             YY             '); 
c = colorbar('Position',[0.906593759638121 0.576208178438662 0.0117370892018781 0.142503097893432]); 
ylabel(c, 'E');
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp
caxis([-1 1])
plotm(study_area(:,2),study_area(:,1),'b','LineWidth',1.5),
subplot(2, 2, 3);
axesm('ortho','origin',[64 -22]), axis off
surfm(LAT, LON, ZZA, topo_int*1e-5-1); 
camlight, material dull, lighting Gouraud, shading interp 
title('(c)           ZZ              '); 
c = colorbar('Position',[0.476857801309071 0.12639405204461 0.0146263836057707 0.166047087980174]); 
ylabel(c, 'E');
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp
caxis([-1 1])
plotm(study_area(:,2),study_area(:,1),'b','LineWidth',1.5),
subplot(2, 2, 4);
axesm('ortho','origin',[64 -22]), axis off
surfm(LAT, LON, XZA, topo_int*1e-5-1); 
camlight, material dull, lighting Gouraud, shading interp
title('(d)          XZ          '); 
c = colorbar('Position',[0.915317329769371 0.122676579925651 0.0126173708920188 0.172242874845105]); 
ylabel(c, 'E');
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp
caxis([-1 1])
plotm(study_area(:,2),study_area(:,1),'b','LineWidth',1.5),
set(gcf,'units','normalized','outerposition',[0.2635 0.1000 0.4365 0.7500])
print('-dpng','-r400','../fig/fig1')
print('-depsc','-r400','../fig/fig1')