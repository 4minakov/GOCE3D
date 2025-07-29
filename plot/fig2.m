%%
%clear all 
%close all

addpath ../data
load GSHHS_i %coastlines

%%
minlat = 55; maxlat = 82; nlat=maxlat-minlat;
minlon = -55; maxlon = 20; nlon=maxlon-minlon;
study_area = [linspace(minlon,maxlon,100)',maxlat*ones(100,1);...
    linspace(maxlon,minlon,100)',minlat*ones(100,1);...
    minlon,maxlat];
%%
% Bathymetry from grid
load topography
% Ice grid
load ice
% Sediment thickness grid
load sediments
% Moho depth (NAGTEC & Szwillus et al. combined)
load Moho_combined
% Age grid
load OceanAge
age1=age;age1(age<0)=max(age(:));
%%
ny0=1081;
nx0=2161;
lontop = linspace(-180,180,nx0);
lattop = linspace(  90,-90,ny0);
[Phi,Lam] = meshgrid(lattop,lontop);
%%
figure(3),clf
%
subplot(2, 2, 1);
ax1 = axesm('ortho','origin',[64 -22]); axis off
surfm(Phi, Lam, topoi'); title('(a)         Topography          '); 
c = colorbar('Position',[0.475641256296906 0.571251548946716 0.0117370892018775 0.151177199504337]); 
ylabel(c, 'm'); 
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp, 
caxis([-4e3 4e3])
plotm(study_area(:,2),study_area(:,1),'r','LineWidth',1.5),
%
subplot(2, 2, 2);
ax2 = axesm('ortho','origin',[64 -22]); axis off
surfm(Phi, Lam, 1e-3*sed_comb'); title('(b)         Sediment thickness          '); 

c = colorbar('Position',[0.906593759638121 0.576208178438662 0.0117370892018781 0.142503097893432]); 
ylabel(c, 'km'); 
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp, 
caxis([1 15])
plotm(study_area(:,2),study_area(:,1),'r','LineWidth',1.5),
%
subplot(2, 2, 3);
ax3 = axesm('ortho','origin',[64 -22]); axis off
surfm(Phi, Lam, 1e-3*M_comb'); title('(c)           Moho depth          '); 
c = colorbar('Position',[0.476857801309071 0.12639405204461 0.0146263836057707 0.166047087980174]);
ylabel(c, 'km'); caxis([10 45])
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp, 
plotm(study_area(:,2),study_area(:,1),'r','LineWidth',1.5),
%
subplot(2, 2, 4);
ax4 = axesm('ortho','origin',[64 -22]); axis off
surfm(Phi, Lam, age1'); 
title('(d)          Oceanic crustal age         '); 
c = colorbar('Position',[0.915317329769371 0.122676579925651 0.0126173708920188 0.172242874845105]); 
ylabel(c, 'Ma'); 
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp, colormap(ax4,flipud(hot))
caxis([0 180])
plotm(study_area(:,2),study_area(:,1),'r','LineWidth',1.5),
%
set(gcf,'units','normalized','outerposition',[0.2635 0.1000 0.4365 0.7500])
%
print('-dpng','-r400','../fig/fig2')
print('-depsc','-r400','../fig/fig2')
