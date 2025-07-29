%% Figure 23. Model variance 

%load figData 
load transects2
load GOCE_NEATLANTIC
load vik
[ny,nx,nr]=size(GOCE_NATL.m);

[Lon3d,Lat3d,r3d]=meshgrid(GOCE_NATL.lon,GOCE_NATL.lat,GOCE_NATL.r);
cm3d = reshape(diag(GOCE_NATL.Cm),ny,nx,nr );

% jj=5;
% 
% cm2d = cm3d(:,:,jj);
% 
% figure,
% ax2=axesm('lambertstd','MapLatLimit',[minlat maxlat-2],'MapLonLimit',[minlon maxlon-8],...
%     'MLineLocation', 10, 'PLineLocation', 5,...
%     'FontSize', 12, 'FontWeight','Bold','LabelFormat','none',...
%     'LabelRotation','on','GLineWidth',.5,'GLineStyle','-','GAltitude',4,...
%     'MLabelParallel',minlat-.01,'PLabelMeridian',minlon-.1);axis off; framem on; gridm off; mlabel on; plabel on;
% surfm(GOCE_NATL.lat,GOCE_NATL.lon,sqrt(cm2d));
% camlight, material dull, lighting Gouraud, shading interp, alpha(1.)
% c=colorbar(ax2,'Position',[0.7 0.65 0.015 0.15]); ylabel(c, '1SD [kg m^{-3}]'); 
% hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp, colormap(ax2,parula(10))
% plotm(transects2(1).lat,transects2(1).lon,'Color',[0.6 0.7 0.6],'LineWidth',2)
% textm(transects2(1).lat(1)-0.1,transects2(1).lon(1)+0.1,'1','Color',[0 0 0],'FontSize',12,'FontWeight','bold')
% plotm(transects2(2).lat,transects2(2).lon,'Color',[0.6 0.7 0.6],'LineWidth',2)
% textm(transects2(2).lat(1)-0.1,transects2(2).lon(1)+0.1,'2','Color',[0 0 0],'FontSize',12,'FontWeight','bold')
% plotm(transects2(3).lat,transects2(3).lon,'Color',[0.6 0.7 0.6],'LineWidth',2)
% textm(transects2(3).lat(1)-0.1,transects2(3).lon(1)+0.1,'3','Color',[0 0 0],'FontSize',12,'FontWeight','bold')
% 
% title([num2str(160),' km'])
% set(gcf,'units','normalized','outerposition',[0.1 0.1 .5 .5])

%%

pr1cm=zeros(nr,length(transects2(1).d1));
pr2cm=zeros(nr,length(transects2(2).d1));
pr3cm=zeros(nr,length(transects2(3).d1));
for i = 1:nr
     pr1cm(i,:) = interp2(GOCE_NATL.lon,GOCE_NATL.lat,cm3d(:,:,i),transects2(1).lon1,transects2(1).lat1,'linear',0);  
     pr2cm(i,:) = interp2(GOCE_NATL.lon,GOCE_NATL.lat,cm3d(:,:,i),transects2(2).lon1,transects2(2).lat1,'linear',0);  
     pr3cm(i,:) = interp2(GOCE_NATL.lon,GOCE_NATL.lat,cm3d(:,:,i),transects2(3).lon1,transects2(3).lat1,'linear',0);  
end
pr1cm = interp1(6371-1e-3*GOCE_NATL.r,sqrt(pr1cm),tomoNA.depth);
pr2cm = interp1(6371-1e-3*GOCE_NATL.r,sqrt(pr2cm),tomoNA.depth);
pr3cm = interp1(6371-1e-3*GOCE_NATL.r,sqrt(pr3cm),tomoNA.depth);
%

figure,
subplot(311)
[d2d,z2d]=meshgrid(transects2(1).d1, tomoNA.depth);
contourf(d2d,z2d, pr1cm), c=colorbar;
hold on, plot(transects2(1).dist,1e-3*transects2(1).LAB,'b--','LineWidth',1)
colormap(flipud(vik)),ylabel(c,'1SD [kg m^{-3}]^2'); title('(a)             '),
ylabel('Depth (km)')
xlim([0 3000]),ylim([0,275]),colormap(gray),caxis([10 20])
tt=title('(a)'); tt.Units='Normalized'; tt.Position(1)=0;
set(gca,'Ydir','reverse')
subplot(312)
[d2d,z2d]=meshgrid(transects2(2).d1, tomoNA.depth);
contourf(d2d,z2d, pr2cm), c=colorbar;
hold on, plot(transects2(1).dist,1e-3*transects2(1).LAB,'b--','LineWidth',1)
ylabel(c,'1SD [kg m^{-3}]'); 
ylabel('Depth (km)'),
xlim([0 3000]),ylim([0,275]),colormap(gray),caxis([10 20])
tt=title('(b)'); tt.Units='Normalized'; tt.Position(1)=0;
set(gca,'Ydir','reverse')
subplot(313)
[d2d,z2d]=meshgrid(transects2(3).d1, tomoNA.depth);
contourf(d2d,z2d, pr3cm), c=colorbar;
hold on, plot(transects2(1).dist,1e-3*transects2(1).LAB,'b--','LineWidth',1)
ylabel(c,'1SD [kg m^{-3}]^2');
ylabel('Depth (km)'),
xlim([0 2500]),ylim([0,275]),colormap(gray),caxis([10 20])
tt=title('(c)'); tt.Units='Normalized'; tt.Position(1)=0;
set(gca,'Ydir','reverse')
set(gcf,'units','normalized','outerposition',[0.1 0.1 .5 .7])
print('-dpng','-r400','../fig/fig23')
print('-depsc','-painters','-r300','../fig/fig23')