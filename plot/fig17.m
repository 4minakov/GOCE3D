%% Figure 17. 3D view models based on inversion of Trr and Txx+Tzz+Tzx
figure,
load real_data_Trr_2_results
load vik
subplot(3,2,1)
T = mLSQR;
plot_spher, colormap(flipud(vik)), %caxis([-10 12])
tt=title('(a)');tt.Units='normalized';
tt.Position(1)=0; tt.HorizontalAlignment='left';
ishore = find([shorelines.Lon] > minlon & [shorelines.Lon] < maxlon & ...
    [shorelines.Lat] > minlat & [shorelines.Lat] < maxlat);
sh_lon = [shorelines.Lon];
sh_lat = [shorelines.Lat];
[sh_x,sh_y,sh_z]=sph2cart(sh_lon(ishore)*pi/180,sh_lat(ishore)*pi/180,sh_lat(ishore)*0+6375e3);
plot3(sh_x,sh_y,sh_z,'.k')
ipb = find(plon > minlon & plon < maxlon & plat > minlat & plat < maxlat);
[pb_x,pb_y,pb_z]=sph2cart(plon(ipb)*pi/180,plat(ipb)*pi/180,plat(ipb)*0+6375e3);
plot3(pb_x,pb_y,pb_z,'.k'), caxis([-20 20]), %c=colorbar('south');
view(90,10),camlight('headlight')
subplot(3,2,3)
T = mSVD;
plot_spher,  %caxis([-10 12])
tt=title('(b)');tt.Units='normalized';
tt.Position(1)=0; tt.HorizontalAlignment='left';
ishore = find([shorelines.Lon] > minlon & [shorelines.Lon] < maxlon & ...
    [shorelines.Lat] > minlat & [shorelines.Lat] < maxlat);
sh_lon = [shorelines.Lon];
sh_lat = [shorelines.Lat];
[sh_x,sh_y,sh_z]=sph2cart(sh_lon(ishore)*pi/180,sh_lat(ishore)*pi/180,sh_lat(ishore)*0+6375e3);
plot3(sh_x,sh_y,sh_z,'.k')
ipb = find(plon > minlon & plon < maxlon & plat > minlat & plat < maxlat);
[pb_x,pb_y,pb_z]=sph2cart(plon(ipb)*pi/180,plat(ipb)*pi/180,plat(ipb)*0+6375e3);
plot3(pb_x,pb_y,pb_z,'.k'), caxis([-20 20]), %c=colorbar('south');
view(90,10),camlight('headlight')
subplot(3,2,5)
T = mPoint;
plot_spher,%caxis([-10 12])
tt=title('(c)');tt.Units='normalized';
tt.Position(1)=0; tt.HorizontalAlignment='left';
ishore = find([shorelines.Lon] > minlon & [shorelines.Lon] < maxlon & ...
    [shorelines.Lat] > minlat & [shorelines.Lat] < maxlat);
sh_lon = [shorelines.Lon];
sh_lat = [shorelines.Lat];
[sh_x,sh_y,sh_z]=sph2cart(sh_lon(ishore)*pi/180,sh_lat(ishore)*pi/180,sh_lat(ishore)*0+6375e3);
plot3(sh_x,sh_y,sh_z,'.k')
ipb = find(plon > minlon & plon < maxlon & plat > minlat & plat < maxlat);
[pb_x,pb_y,pb_z]=sph2cart(plon(ipb)*pi/180,plat(ipb)*pi/180,plat(ipb)*0+6375e3);
plot3(pb_x,pb_y,pb_z,'.k'), caxis([-20 20]), %c=colorbar('south');
view(90,10),camlight('headlight')
%
%
%
load real_data_TzzTxxTzx_2_results
%
subplot(3,2,2)
T = mLSQR;
plot_spher, %caxis([-10 12])
tt=title('(d)');tt.Units='normalized';
tt.Position(1)=0; tt.HorizontalAlignment='left';
ishore = find([shorelines.Lon] > minlon & [shorelines.Lon] < maxlon & ...
    [shorelines.Lat] > minlat & [shorelines.Lat] < maxlat);
sh_lon = [shorelines.Lon];
sh_lat = [shorelines.Lat];
[sh_x,sh_y,sh_z]=sph2cart(sh_lon(ishore)*pi/180,sh_lat(ishore)*pi/180,sh_lat(ishore)*0+6375e3);
plot3(sh_x,sh_y,sh_z,'.k')
ipb = find(plon > minlon & plon < maxlon & plat > minlat & plat < maxlat);
[pb_x,pb_y,pb_z]=sph2cart(plon(ipb)*pi/180,plat(ipb)*pi/180,plat(ipb)*0+6375e3);
plot3(pb_x,pb_y,pb_z,'.k'), caxis([-20 20]), %c=colorbar('south');
view(90,10),camlight('headlight')
subplot(3,2,4)
T = mSVD;
plot_spher,  %caxis([-10 12])
tt=title('(e)');tt.Units='normalized';
tt.Position(1)=0; tt.HorizontalAlignment='left';
ishore = find([shorelines.Lon] > minlon & [shorelines.Lon] < maxlon & ...
    [shorelines.Lat] > minlat & [shorelines.Lat] < maxlat);
sh_lon = [shorelines.Lon];
sh_lat = [shorelines.Lat];
[sh_x,sh_y,sh_z]=sph2cart(sh_lon(ishore)*pi/180,sh_lat(ishore)*pi/180,sh_lat(ishore)*0+6375e3);
plot3(sh_x,sh_y,sh_z,'.k')
ipb = find(plon > minlon & plon < maxlon & plat > minlat & plat < maxlat);
[pb_x,pb_y,pb_z]=sph2cart(plon(ipb)*pi/180,plat(ipb)*pi/180,plat(ipb)*0+6375e3);
plot3(pb_x,pb_y,pb_z,'.k'), caxis([-20 20]), %c=colorbar('south');
view(90,10),camlight('headlight')
subplot(3,2,6)
T = mPoint;
plot_spher, %caxis([-10 12])
tt=title('(f)');tt.Units='normalized';
tt.Position(1)=0; tt.HorizontalAlignment='left';
ishore = find([shorelines.Lon] > minlon & [shorelines.Lon] < maxlon & ...
    [shorelines.Lat] > minlat & [shorelines.Lat] < maxlat);
sh_lon = [shorelines.Lon];
sh_lat = [shorelines.Lat];
[sh_x,sh_y,sh_z]=sph2cart(sh_lon(ishore)*pi/180,sh_lat(ishore)*pi/180,sh_lat(ishore)*0+6375e3);
plot3(sh_x,sh_y,sh_z,'.k')
ipb = find(plon > minlon & plon < maxlon & plat > minlat & plat < maxlat);
[pb_x,pb_y,pb_z]=sph2cart(plon(ipb)*pi/180,plat(ipb)*pi/180,plat(ipb)*0+6375e3);
plot3(pb_x,pb_y,pb_z,'.k'), caxis([-20 20]), %c=colorbar('south');
view(90,10),camlight('headlight')
%
%
c=colorbar('south');c.Position=[0.4 0.07 0.25 0.014];
c.FontSize=10; ylabel(c,'kg m$^{-3}$','Interpreter','latex','FontSize',12)
%
set(gcf,'Units','normalized','OuterPosition',[0.1 0.1 .5 0.75])
print('-dpng','-r400','../fig/fig17')
print('-depsc','-r400','../fig/fig17')