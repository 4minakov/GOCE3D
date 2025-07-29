% figure 19
% Density and dVs variation extracted along Profile 1
% extract model values along transects
load transects2
[nxs,nys,nzs]=size(tomoNA.vsh);
pr1vs=zeros(nzs,length(transects2(1).d1));
for iz = 1:nzs
 pr1vs(iz,:) = interp2(tomoNA.longitude, tomoNA.latitude, tomoNA.vsh(:,:,iz)',transects2(1).lon1,transects2(1).lat1,'linear',0);
end
pr1vs0 = repmat(squeeze(mean(mean(tomoNA.vsh,1),2)),1,length(transects2(1).d1));
dvs1 = (pr1vs-pr1vs0)./pr1vs0 * 100;

% mantle density from rifting model
z_depth = zeros(Model.number_of_layers,1);
dens_td = zeros(180,360,Model.number_of_layers);
pr1rht=zeros(Model.number_of_layers,length(transects2(1).d1));
lon1 = transects2(1).lon1; lon1(lon1<0)=lon1(lon1<0)+360;
for i=1:Model.number_of_layers
     Lupper = Model.(['l',num2str(i)]);
     z_depth(i) = -1e-3*mean(mean(Lupper.bound));
     dens_td(:,:,i) = Lupper.dens + zeros(180,360);
     pr1rht(i,:) = interp2(Model.Lon,Model.Lat,dens_td(:,:,i),lon1,transects2(1).lat1);
 end
pr1rht = interp1(z_depth,pr1rht,tomoNA.depth);
pr1rht0 = repmat(squeeze(mean(mean(dens_td(140:end,310:end,:),1),2)),1,length(transects2(1).d1));
pr1rht0 = interp1(z_depth,pr1rht0,tomoNA.depth);
drhot1 = (pr1rht-pr1rht0); drhot1(isnan(drhot1))=0;

% density model 3D gravity inversion 
[ny,nx,nr]=size(rhoi.m3d);
pr1rh=zeros(nr,length(transects2(1).d1));
for i = 1:nr
     pr1rh(i,:) = interp2(rhoi.Lon_reg,rhoi.Lat_reg,rhoi.m3d(:,:,i),transects2(1).lon1,transects2(1).lat1,'linear',0);  
end
pr1rh = interp1(6371-1e-3*rhoi.r,pr1rh,tomoNA.depth);
dpr1rh = pr1rh./pr1rht0 * 100;

%%
figure,
subplot(411), plot(transects2(1).dist, transects2(1).dg), axis tight
ylabel('mGal'),%xlim([-51 -4])
xlim([0 3000])
tt=title('(a)'); tt.Units='Normalized'; tt.Position(1)=0; legend('Free-air gravity','Location','northeast')
set(gca,'Position',[0.15 0.85 .7 .1])
%
subplot(412), plot(transects2(1).dist, transects2(1).moho), hold on
plot(transects2(1).dist, transects2(1).basement), hold on
plot(transects2(1).dist, transects2(1).topo), hold on
xlim([0 3000])
legend('MOHO','BSM','TOPO','Location','west')
set(gca,'Position',[0.15 .60 .7 .2])
set(gca,'Ydir','reverse'),ylabel('Depth (km)'), 
tt=title('(b)'); tt.Units='Normalized'; tt.Position(1)=0;
%
subplot(413), 
imagesc(transects2(1).d1, tomoNA.depth, pr1rh+drhot1), c=colorbar;
hold on, plot(transects2(1).dist,1e-3*transects2(1).LAB,'k--','LineWidth',1)
colormap(flipud(rwbcmap)),ylabel(c,'[kg m^{-3}]'); title('(c)             '),
ylabel('Depth (km)'),%xlim([-51 -4])
xlim([0 3000]),ylim([50,300])
tt=title('(c)'); tt.Units='Normalized'; tt.Position(1)=0;
set(gca,'Position',[0.15 .34 .7 .2]) 
caxis([-50 50])
%
subplot(414), 
imagesc(transects2(1).d1, tomoNA.depth,dvs1 ),  
hold on, plot(transects2(1).dist,1e-3*transects2(1).LAB,'k--','LineWidth',1)
c=colorbar;
colormap(flipud(rwbcmap)),ylabel(c,'[%]'); title('(d)                   '),
xlabel('Distance (km)'),ylabel('Depth (km)'), 
shading interp,
xlim([0 3000]),ylim([50 300])
tt=title('(d)'); tt.Units='Normalized'; tt.Position(1)=0;
set(gca,'Position',[0.15 .07 .7 .2]),caxis([-10 10]) 
%
set(gcf,'units','normalized','outerposition',[0.1 0.1 .5 .5])
print('-dpng','-r400','../fig/fig19')
print('-depsc','-painters','-r300','../fig/fig19')
