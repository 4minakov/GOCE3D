%% data covariance from moho variance
COT=[]
modfolder = 'C:\Users\alexamin\Dropbox (UiO)\ESA2\models\';% path to folder with models
addpath \\kant\geo-ceed-u1\alexamin\3DEarth\
addpath \\kant\geo-ceed-u1\alexamin\gravity
load NA_DATASET
load ResidualAnomalyNA
load ResidualTrrNA
COT=[]; %plotting flag
% parameters of polar stereographic projection
phi0 =  65;%true latitude
lam0 = -70;%zero longitude
dxm = 30e3;
dym = 30e3;
xmmin = -0;
xmmax =  4000e3;
ymmin = -1000e3;
ymmax =  3500e3;
% Moho based on refraction data in CRUST1.0 from W.Swillius
%m_globe = load([modfolder,'crust\Moho-interp-06-04-2017.txt']);
m_globe = load('C:\Users\alexamin\Dropbox (UiO)\ESA2\data\crust\szwillus\jgrb53251-sup-0003-data_set_si-s03.txt');
m_km = reshape(m_globe(:,3),360,180);
lon_m = reshape(m_globe(:,1),360,180);
lat_m = reshape(m_globe(:,2),360,180);
err_m = reshape(m_globe(:,4),360,180);
crust1_pt = load([modfolder,'crust\gsc-only-locations.csv']);
[Xpt,Ypt]= pstereo(crust1_pt(:,2),crust1_pt(:,1),phi0,lam0);
[X_crust1, Y_crust1] = pstereo(lat_m,lon_m,phi0,lam0);
%
FCR = scatteredInterpolant(X_crust1(:),Y_crust1(:),m_km(:)); 
FCR.Method = 'natural';
FCR.ExtrapolationMethod='none';
moho = FCR(Xmesh*1e-3,Ymesh*1e-3);
%
moho_pt = FCR(Xpt, Ypt);
ii = find(crust1_pt(:,1)>-60&crust1_pt(:,1)<30&crust1_pt(:,2)>50&crust1_pt(:,2)<87);
ii1 = find((lon_m(:)<-60&lon_m(:)>-80&lat_m(:)<87&lat_m(:)>83) | ...
    (lon_m(:)>35&lon_m(:)<45&lat_m(:)>40&lat_m(:)<55)  );


Fpt = scatteredInterpolant(Xpt(ii),Ypt(ii),moho_pt(ii));
Fpt.Method = 'natural';
Fpt.ExtrapolationMethod='none';
moho_pt_grd = Fpt(Xmesh*1e-3,Ymesh*1e-3);
%
FCR1 = scatteredInterpolant(X_crust1(:),Y_crust1(:),err_m(:)); 
FCR1.Method = 'natural';
FCR1.ExtrapolationMethod='none';
Err_moho = FCR1(Xmesh*1e-3,Ymesh*1e-3);
%
%fid = fopen('CRUST1_Moho_NA_pt.txt','w');
%fprintf(fid,'%9.4f,%9.4f,%9.1f\n',[crust1_pt(ii,1)';crust1_pt(ii,2)';-1000*moho_pt(ii)']);
%fclose(fid);
%% load kriging results usiing Wolfgang's program
mkrig = load('test2.out');
mkrig_var = load('test2_var.out');
[Xkrig,Ykrig]= pstereo(mkrig(:,2),mkrig(:,1),phi0,lam0);
FKRIG = scatteredInterpolant(Xkrig,Ykrig,mkrig(:,3)); 
FKRIG.Method = 'natural';
FKRIG.ExtrapolationMethod='none';
mkrig_grd = FKRIG(Xmesh*1e-3,Ymesh*1e-3);
FKRIG.Values=mkrig_var(:,3);
mkrigvar_grd = FKRIG(Xmesh*1e-3,Ymesh*1e-3);

%% Nagtec profiles

%%
%load oleron
%load bamako 
addpath ./data
fid = 'C:\Users\alexamin\Dropbox (UiO)\ESA2\data\crust\nagtec\Profiles\NAG-TEC_refraction.xy';

pr = xygmt_load(fid);
xy1 = [];
ellipsoid = wgs84Ellipsoid;
for i = 1:200
    xyn = pr.(['l',num2str(i)]);
    for j=1:size(xyn,1)-1
        dst = distance(xyn(j,2),xyn(j,1),xyn(j+1,2),xyn(j+1,1),ellipsoid,'degree');
        npts = round(dst/30e3)+1;
        [lattrk,lontrk] = track(xyn(:,2),xyn(:,1),ellipsoid,'degree',npts);
        xy1 = [xy1; [lattrk, lontrk] ];
    end
%     figure(36),
%     [prx,pry]=pstereo(xyn(:,2),xyn(:,1),phi0,lam0);
%     text(prx+10,pry+10,num2str(i)), drawnow
end
[X_nagtec,Y_nagtec]= pstereo(xy1(:,1),xy1(:,2),phi0,lam0);
% 
% 
% 
% FNAG = scatteredInterpolant(Xmesh(:)*1e-3,Ymesh(:)*1e-3,NA_DATASET.Mo(:)); 
% FNAG.Method = 'natural';
% FNAG.ExtrapolationMethod='none';
% mNAG = FNAG(X_nagtec,Y_nagtec);
%
%
modfolder = 'C:\Users\alexamin\Dropbox (UiO)\ESA2\models\';
src1 = [modfolder,'crust\nagtec\depth_to_moho_from_seismic_refraction_geo.nc']
finfo = ncinfo(src1)
X1 = ncread(src1,'x');
Y1 = ncread(src1,'y');
M_nag  = ncread(src1,'z'); 
%M_nag(isnan(M_nag))=0;
%M_nag=M_nag';
[Lon_nag,Lat_nag] = meshgrid(X1,Y1);
%Lon_nag(isnan(M_nag))=[];
%Lat_nag(isnan(M_nag))=[];
%M_nag(isnan(M_nag))=[];
[X_nag, Y_nag] = pstereo(Lat_nag,Lon_nag,phi0,lam0);
%F1 = scatteredInterpolant(X_nag(:), Y_nag(:), M_nag(:));
F1 = griddedInterpolant({X1, Y1}, M_nag);
F1.Method = 'linear';
F1.ExtrapolationMethod='none';
mNAG = F1(xy1(:,2),xy1(:,1));
mNAG_GlobMesh = F1(lon_m,lat_m);
%%
figure, scatter(X_nagtec,Y_nagtec,5,mNAG), hold on
hold on, scatter(crust1_pt(ii,1),crust1_pt(ii,2),5,moho_pt(ii))
%ii2 = find( X_crust1<xmmin|X_crust1>xmmax|Y_crust1<ymmin|Y_crust1>ymmax) ;
xcr = X_crust1(isnan(mNAG_GlobMesh)); 
ycr = Y_crust1(isnan(mNAG_GlobMesh)); 
zcr = m_km(isnan(mNAG_GlobMesh));
loncr = lon_m(isnan(mNAG_GlobMesh));
latcr = lat_m(isnan(mNAG_GlobMesh));
ixmask = find(xcr>xmmin*1e-3&xcr<xmmax*1e-3&ycr>ymmin*1e-3&ycr<ymmax*1e-3);
scatter(xcr(ixmask(1:5:end)),ycr(ixmask(1:5:end)),5,zcr(ixmask(1:5:end)))

%plot_basemap
%% combined
% mm_out = ...
% [crust1_pt(ii,1)', xy1(~isnan(mNAG),2)';...
%  crust1_pt(ii,2)', xy1(~isnan(mNAG),1)';...
% -1000*moho_pt(ii)', -1000*mNAG(~isnan(mNAG))' ];
mm_out = ...
[xy1(:,2),xy1(:,1),-1e3*mNAG;...
 crust1_pt(ii,1),crust1_pt(ii,2),-1e3*moho_pt(ii);...
 loncr(ixmask(1:5:end)),latcr(ixmask(1:5:end)),-1e3*zcr(ixmask(1:5:end))];
mm_out(isnan(mm_out(:,3)),:)=[];
mm_out=mm_out';
%fid = fopen('CRUST1_NAGTEC_Moho_NA_pt2.txt','w');
%fprintf(fid,'%9.4f,%9.4f,%9.1f\n',mm_out);
%fclose(fid);
%%
load haxby_GMT.mat
figure(30),clf
surf(Xmesh*1e-3,Ymesh*1e-3,NA_DATASET.Topo*3e-2-100,moho), 
camlight, material dull, lighting Gouraud, shading interp
colormap(haxby_GMT), caxis([1 40])
view(0,90), axis equal tight, camlight, shading interp,hold on
cb = colorbar('Location','southoutside'); cb.Position=[0.4 0.05 0.2 0.05];
cb.Label.String = 'km'; cb.Label.Position = [45 1];
hold on, title('Global grid')
plot(Xpt(ii),Ypt(ii),'o','MarkerFaceColor','k','MarkerSize',3)
plot_basemap
%
figure(31),clf
surf(Xmesh*1e-3,Ymesh*1e-3,NA_DATASET.Topo*3e-2-100,moho_pt_grd), 
camlight, material dull, lighting Gouraud, shading interp
colormap(haxby_GMT), caxis([1 40])
view(0,90), axis equal tight, camlight, shading interp,hold on
cb = colorbar('Location','southoutside'); cb.Position=[0.4 0.05 0.2 0.05];
cb.Label.String = 'km'; cb.Label.Position = [45 1];
hold on, title('Interpolation of point data')
plot(Xpt(ii),Ypt(ii),'o','MarkerFaceColor','k','MarkerSize',3)
plot_basemap

figure(32),clf
surf(Xmesh*1e-3,Ymesh*1e-3,NA_DATASET.Topo*3e-2-100,moho_pt_grd-moho), 
camlight, material dull, lighting Gouraud, shading interp
colormap(haxby_GMT), %caxis([1 40])
view(0,90), axis equal tight, camlight, shading interp,hold on
cb = colorbar('Location','southoutside'); cb.Position=[0.4 0.05 0.2 0.05];
cb.Label.String = 'km'; cb.Label.Position = [45 1];
hold on, title('Difference')
plot(Xpt(ii),Ypt(ii),'o','MarkerFaceColor','k','MarkerSize',3)
plot_basemap
%%
figure(331),clf
surf(Xmesh*1e-3,Ymesh*1e-3,NA_DATASET.Topo*3e-2-100,mkrig_grd), 
camlight, material dull, lighting Gouraud, shading interp
colormap(colbrew(4,12)), caxis([1 40])
view(0,90), axis equal tight, camlight, shading interp,hold on
cb = colorbar('Location','southoutside'); cb.Position=[0.4 0.05 0.2 0.05];
cb.Label.String = 'km'; cb.Label.Position = [45 1];
cb.FontSize=12;
hold on, %title('Kriging Moho')
plot(Xpt(ii),Ypt(ii),'o','MarkerFaceColor','k','MarkerEdgeColor','g','MarkerSize',3)
plot(X_nagtec,Y_nagtec,'o','MarkerFaceColor','k','MarkerEdgeColor','g','MarkerSize',3)
plot_basemap
set(gcf,'units','normalized')
set(gcf,'Position',[0.0635 0.2808 0.4594 0.6483]);
print('-depsc','-painters','-r400','Fig_Moho_krig')
print('-dpng','-r400','Fig_Moho_krig')
%%
figure(341),clf
surf(Xmesh*1e-3,Ymesh*1e-3,NA_DATASET.Topo*3e-2-100,sqrt(mkrigvar_grd)), 
camlight, material dull, lighting Gouraud, shading interp, alpha(0.9)
colormap(gray(12)), %([1.5 4])
view(0,90), axis equal tight, camlight, shading interp,hold on
cb = colorbar('Location','southoutside'); cb.Position=[0.4 0.05 0.2 0.05];
cb.Label.String = 'km'; cb.Label.Position = [6.7 0.8]; %cb.Label.FontSize=12;
cb.FontSize=12; cb.Position = [0.672094035698744 0.145120571751759 0.13974386881439 0.0310324438274476];
hold on, %title('Kriging Moho Std')
plot(Xpt(ii),Ypt(ii),'o','MarkerFaceColor','k','MarkerEdgeColor','g','MarkerSize',3)
plot(X_nagtec,Y_nagtec,'o','MarkerFaceColor','k','MarkerEdgeColor','g','MarkerSize',3)
plot_basemap
set(gcf,'units','normalized')
set(gcf,'Position',[0.0635 0.2808 0.4594 0.6483]);
%print('-depsc','-painters','-r400','Fig_Moho_krigstd')
print('-dpng','-r400','Fig_Moho_krigstd3')
%%
figure(35),clf
surf(Xmesh*1e-3,Ymesh*1e-3,NA_DATASET.Topo*3e-2-100,Err_moho), 
camlight, material dull, lighting Gouraud, shading interp
colormap(haxby_GMT), %caxis([1 40])
view(0,90), axis equal tight, camlight, shading interp,hold on
cb = colorbar('Location','southoutside'); cb.Position=[0.4 0.05 0.2 0.05];
cb.Label.String = 'km'; cb.Label.Position = [45 1];
hold on, title('Global grid Std')
plot(Xpt(ii),Ypt(ii),'o','MarkerFaceColor','k','MarkerSize',3)
plot_basemap

figure(36),clf
surf(Xmesh*1e-3,Ymesh*1e-3,NA_DATASET.Topo*3e-2-100,NA_DATASET.Mo), 
camlight, material dull, lighting Gouraud, shading interp
colormap(haxby_GMT), caxis([1 40])
view(0,90), axis equal tight, camlight, shading interp,hold on
cb = colorbar('Location','southoutside'); cb.Position=[0.4 0.05 0.2 0.05];
cb.Label.String = 'km'; cb.Label.Position = [45 1];
hold on, title('Moho NAGTEC')
plot(Xpt(ii),Ypt(ii),'o','MarkerFaceColor','k','MarkerSize',3)
plot(X_nagtec,Y_nagtec,'o','MarkerFaceColor','r','MarkerSize',3)
plot_basemap


figure(37),clf
surf(Xmesh*1e-3,Ymesh*1e-3,NA_DATASET.Topo*3e-2-100,mkrig_grd-NA_DATASET.Mo), 
camlight, material dull, lighting Gouraud, shading interp
colormap(haxby_GMT), caxis([-10 10])
view(0,90), axis equal tight, camlight, shading interp,hold on
cb = colorbar('Location','southoutside'); cb.Position=[0.4 0.05 0.2 0.05];
cb.Label.String = 'km'; cb.Label.Position = [45 1];
hold on, title('Kriging Moho Std')
plot(Xpt(ii),Ypt(ii),'o','MarkerFaceColor','k','MarkerEdgeColor','b','MarkerSize',3)
plot(X_nagtec,Y_nagtec,'o','MarkerFaceColor','k','MarkerEdgeColor','b','MarkerSize',3)
plot_basemap

%%
lon_krig2d = reshape(mkrig(:,1),85,66);
lat_krig2d = reshape(mkrig(:,2),85,66);
m_krig2d = reshape(mkrig(:,3),85,66);
[xkr2d,ykr2d]=pstereo(lat_krig2d,lon_krig2d,phi0,lam0);
figure(38),clf
surf(xkr2d,ykr2d,m_krig2d), 
camlight, material dull, lighting Gouraud, shading interp
colormap(haxby_GMT), caxis([10 40])
view(0,90), axis equal tight, camlight, shading interp,hold on
cb = colorbar('Location','southoutside'); cb.Position=[0.4 0.05 0.2 0.05];
cb.Label.String = 'km'; cb.Label.Position = [45 1];
hold on, title('Kriging Moho Std')
plot(Xpt(ii),Ypt(ii),'o','MarkerFaceColor','k','MarkerEdgeColor','b','MarkerSize',3)
plot(X_nagtec,Y_nagtec,'o','MarkerFaceColor','k','MarkerEdgeColor','b','MarkerSize',3)
plot_basemap

figure(39),clf
surf(X_nag,Y_nag,M_nag'), 
%camlight, material dull, lighting Gouraud, shading interp
colormap(haxby_GMT), caxis([10 40])
view(0,90), axis equal tight, %camlight, 
shading interp,hold on
cb = colorbar('Location','southoutside'); cb.Position=[0.4 0.05 0.2 0.05];
cb.Label.String = 'km'; cb.Label.Position = [45 1];
hold on, title('Kriging Moho Std')
plot(Xpt(ii),Ypt(ii),'o','MarkerFaceColor','k','MarkerEdgeColor','b','MarkerSize',3)
plot(X_nagtec,Y_nagtec,'o','MarkerFaceColor','k','MarkerEdgeColor','b','MarkerSize',3)
plot_basemap

%set(gcf,'units','normalized','outerposition',[0.1 0.1 .3 .5])
%print('-dpng','-r300','-opengl','MohoCR1')