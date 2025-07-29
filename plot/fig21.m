% figure 21
% Density and dVs variation extracted along Profile 3
% extract model values along transects
load transects2

pr3vs=zeros(nzs,length(transects2(3).d1));
for iz = 1:nzs
 pr3vs(iz,:) = interp2(tomoNA.longitude, tomoNA.latitude, tomoNA.vsh(:,:,iz)',transects2(3).lon1,transects2(3).lat1,'linear',0);
end
pr3vs0 = repmat(squeeze(mean(mean(tomoNA.vsh,1),2)),1,length(transects2(3).d1));
dvs3 = (pr3vs-pr3vs0)./pr3vs0 * 100;

% mantle density from rifting model
z_depth = zeros(Model.number_of_layers,1);
dens_td = zeros(180,360,Model.number_of_layers);
pr3rht=zeros(Model.number_of_layers,length(transects2(3).d1));
lon1 = transects2(3).lon1; lon1(lon1<0)=lon1(lon1<0)+360;
for i=1:Model.number_of_layers
     Lupper = Model.(['l',num2str(i)]);
     z_depth(i) = -1e-3*mean(mean(Lupper.bound));
     dens_td(:,:,i) = Lupper.dens + zeros(180,360);
     pr3rht(i,:) = interp2(Model.Lon,Model.Lat,dens_td(:,:,i),lon1,transects2(3).lat1);
 end
pr3rht = interp1(z_depth,pr3rht,tomoNA.depth);
pr3rht0 = repmat(squeeze(mean(mean(dens_td(140:end,310:end,:),1),2)),1,length(transects2(3).d1));
pr3rht0 = interp1(z_depth,pr3rht0,tomoNA.depth);
drhot3 = (pr3rht-pr3rht0); drhot3(isnan(drhot3))=0;

% density model 3D gravity inversion 
[ny,nx,nr]=size(rhoi.m3d);
pr3rh=zeros(nr,length(transects2(3).d1));
for i = 1:nr
     pr3rh(i,:) = interp2(rhoi.Lon_reg,rhoi.Lat_reg,rhoi.m3d(:,:,i),transects2(3).lon1,transects2(3).lat1,'linear',0);  
end
pr3rh = interp1(6371-1e-3*rhoi.r,pr3rh,tomoNA.depth);
dpr3rh = pr3rh./pr3rht0 * 100;
%%
figure,
subplot(411), plot(transects2(3).dist, transects2(3).dg), axis tight
ylabel('mGal'),%xlim([-51 -4])
xlim([0 3000])
tt=title('(a)'); tt.Units='Normalized'; tt.Position(1)=0; legend('Free-air gravity','Location','northeast')
set(gca,'Position',[0.15 0.85 .7 .1])
%
subplot(412), plot(transects2(3).dist, transects2(3).moho), hold on
plot(transects2(3).dist, transects2(3).basement), hold on
plot(transects2(3).dist, transects2(3).topo), hold on
xlim([0 3000])
legend('MOHO','BSM','TOPO','Location','west')
set(gca,'Position',[0.15 .60 .7 .2])
set(gca,'Ydir','reverse'),ylabel('Depth (km)'), 
tt=title('(b)'); tt.Units='Normalized'; tt.Position(1)=0;
%
subplot(413), 
imagesc(transects2(3).d1, tomoNA.depth, pr3rh+drhot3), c=colorbar;
hold on, plot(transects2(3).dist,1e-3*transects2(3).LAB,'k--','LineWidth',1)
colormap(flipud(rwbcmap)),ylabel(c,'[kg m^{-3}]'); title('(c)             '),
ylabel('Depth (km)'),%xlim([-51 -4])
xlim([0 3000]),ylim([50,300])
tt=title('(c)'); tt.Units='Normalized'; tt.Position(1)=0;
set(gca,'Position',[0.15 .34 .7 .2]) 
caxis([-50 50])
%
subplot(414), 
imagesc(transects2(3).d1, tomoNA.depth,dvs3 ),  
hold on, plot(transects2(3).dist,1e-3*transects2(3).LAB,'k--','LineWidth',1)
c=colorbar;
colormap(flipud(rwbcmap)),ylabel(c,'[%]'); title('(d)                   '),
xlabel('Distance (km)'),ylabel('Depth (km)'), 
shading interp,
xlim([0 3000]),ylim([50 300])
tt=title('(d)'); tt.Units='Normalized'; tt.Position(1)=0;
set(gca,'Position',[0.15 .07 .7 .2]),caxis([-10 10]) 
%
set(gcf,'units','normalized','outerposition',[0.1 0.1 .5 .5])
print('-dpng','-r400','../fig/fig21')
print('-depsc','-painters','-r300','../fig/fig21')
