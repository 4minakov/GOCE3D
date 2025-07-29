%clear all 
%close all
%
if ~isfolder('../tools/SHbundle-master')
    disp('Downloading SHbundle m-files ..')
    websave('../tools/SHbundle-master.zip','https://www.gis.uni-stuttgart.de/dokumente/SHbundle-master.zip')
    unzip('../tools/SHbundle-master.zip','../tools/SHbundle-master')
    system('del ..\tools\*.zip')
end
if ~isfolder('../tools/uberall-master')
    disp('Downloading SHbundle m-files ..')
    websave('../tools/uberall-master.zip','https://www.gis.uni-stuttgart.de/dokumente/uberall-master.zip')
    unzip('../tools/uberall-master.zip','../tools/uberall-master')
    system('del ..\tools\*.zip')
end
if ~isfolder('../tools/ScientificColourMaps7')
    disp('Downloading colormaps ..')
    websave('../tools/ScientificColourMaps7.zip','https://zenodo.org/record/5501399/files/ScientificColourMaps7.zip')
    unzip('../tools/ScientificColourMaps7.zip','../tools/ScientificColourMaps7')
    system('del ..\tools\*.zip')
end
if ~isfile('../data/XGM2016.gfc')
    disp('Downloading XGM gravity model ..')
    websave('../data/XGM2016.zip','https://datapub.gfz-potsdam.de/download/10.5880.ICGEM.2017.003/XGM2016.zip')
    unzip('../data/XGM2016.zip','../data')
    system('del ..\data\*.zip')
end
%%
addpath(genpath('../data'))
addpath(genpath('../tools'))
addpath('C:\Program Files\gmt6\bin')
load GSHHS_i
load XGM
load topography
load roma
%% plate boundaries
% filename = 'All_boundaries.txt';
% [plat,plon] = importPlates(filename);
%% XGM2019 from coefficients
gfc_file = 'XGM2016.gfc';
[gsm, lmax, lmin, info] = parse_icgem(gfc_file, 'max_lm', 180); % read model
%%
%
%lmax = XGM.lmax;
%
[field, lmax] = clm2sc(gsm, 'max_lm', lmax); % convert format of coeffs to /S|C\
field1=cs2sc(field); field1(1:3,1:3)=0;
[V_pot, theRAD, lamRAD] = gshs_(field1, 'quant', 'potential', 'grid', 'block', 'gridsize', 2*lmax, 'height', 0, 'sub_wgs84', true);
V_dg = gshs_(field1, 'quant', 'dg', 'grid', 'block', 'gridsize', 2*lmax, 'height', 0, 'sub_wgs84', 1);
V_dg_220 = gshs_(field1, 'quant', 'dg', 'grid', 'block', 'gridsize', 2*lmax, 'height', 225e3, 'sub_wgs84', 1);
V_rr = gshs_(field1, 'quant', 'trr', 'grid', 'block', 'gridsize', 2*lmax, 'height', 0, 'sub_wgs84', 1);
V_rr_220 = gshs_(field1, 'quant', 'trr', 'grid', 'block', 'gridsize', 2*lmax, 'height', 225e3, 'sub_wgs84', 1);
%%
minlat = 55; maxlat = 82; nlat=maxlat-minlat;
minlon = -55; maxlon = 20; nlon=maxlon-minlon;
study_area = [linspace(minlon,maxlon,100)',maxlat*ones(100,1);...
    linspace(maxlon,minlon,100)',minlat*ones(100,1);...
    minlon,maxlat];
%%
[plat,plon] = importPlates('All_boundaries.txt');
ny0=1081; nx0=2161;
lontop = linspace(-180,180,nx0); lattop = linspace(  90,-90,ny0);
[Phi,Lam] = meshgrid(lattop,lontop);
[Phi2d,Lam2d]=meshgrid(90-theRAD*180/pi, [lamRAD,lamRAD(1)]*180/pi);
Lam2d(Lam2d>180)=Lam2d(Lam2d>180)-360;
topo_int = interp2(Phi,Lam,topoi',Phi2d, Lam2d);
%
figure(1),clf
subplot(2, 2, 1);
axesm('ortho','origin',[64 -22]), axis off
surfm(90-theRAD*180/pi, [lamRAD,lamRAD(1)]*180/pi, [V_pot,V_pot(:,1)],1e-5*topo_int'-1); 
camlight, material dull, lighting Gouraud, shading interp
title('(a)                              '); 
c = colorbar('Position',[0.475641256296906 0.571251548946716 0.0117370892018775 0.151177199504337]);
ylabel(c, '[m^2 s^{-2}]'); 
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp, 

addpath(genpath('C:\Users\alexamin\Dropbox (UiO)\Random_Fields_NS\spherical\SHbundle-master'))
colormap(flipud(roma))
colormap(colbrew(1))
plotm(study_area(:,2),study_area(:,1),'b','LineWidth',1.5),
plotm(plat,plon,'k'),
subplot(2, 2, 2);
axesm('ortho','origin',[64 -22]), axis off
surfm(90-theRAD*180/pi, [lamRAD,lamRAD(1)]*180/pi, [V_dg,V_dg(:,1)],1e-5*topo_int'-1); 
camlight, material dull, lighting Gouraud, shading interp
title('(b)                              '); 
c = colorbar('Position',[0.906593759638121 0.576208178438662 0.0117370892018781 0.142503097893432]); 
ylabel(c, '[10^{-5} m s^{-2}]');caxis([-100 100])
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp
plotm(study_area(:,2),study_area(:,1),'b','LineWidth',1.5),
plotm(plat,plon,'k'),
subplot(2, 2, 3);
axesm('ortho','origin',[64 -22]), axis off
surfm(90-theRAD*180/pi, [lamRAD,lamRAD(1)]*180/pi, [V_rr,V_rr(:,1)],1e-5*topo_int'-1); 
camlight, material dull, lighting Gouraud, shading interp
title('(c)                         '); 
c = colorbar('Position',[0.476857801309071 0.12639405204461 0.0146263836057707 0.166047087980174]); 
ylabel(c, '[10^{-9} s^{-2}]');
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp
caxis([-8 8])
plotm(study_area(:,2),study_area(:,1),'b','LineWidth',1.5),
plotm(plat,plon,'k'),
subplot(2, 2, 4);
axesm('ortho','origin',[64 -22]), axis off
%surfm(90-theRAD*180/pi, [lamRAD,lamRAD(1)]*180/pi, [V_rr_220,V_rr_220(:,1)]); title('2nd radial derivative (220 km)'); c = colorbar; ylabel(c, '[1/s^2]');
surfm(90-theRAD*180/pi, [lamRAD,lamRAD(1)]*180/pi, [V_dg_220,V_dg_220(:,1)] ,1e-5*topo_int'-1); 
camlight, material dull, lighting Gouraud, shading interp
title('(d)                             '); 
c = colorbar('Position',[0.915317329769371 0.122676579925651 0.0126173708920188 0.172242874845105]); 
ylabel(c, '[10^{-5} m s^{-2}]');%caxis([-10 30])
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp
%caxis([-1 1])
plotm(study_area(:,2),study_area(:,1),'b','LineWidth',1.5),
plotm(plat,plon,'k'),

set(gcf,'units','normalized','outerposition',[0.2635 0.1000 0.4365 0.7500])
print('-dpng','-r400','../fig/fig3')
print('-depsc','-r400','../fig/fig3')

gmt('psconvert','../fig/fig3.eps -Tf -P -A ')
open('../fig/fig3.pdf')