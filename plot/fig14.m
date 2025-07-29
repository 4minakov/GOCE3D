% figure 14
% Inversion misfit T_rr (predicted minus observed)
figure,
ax2=axesm('lambertstd','MapLatLimit',[minlat maxlat-2],'MapLonLimit',[minlon maxlon-8],...
    'MLineLocation', 10, 'PLineLocation', 5,...
    'FontSize', 12, 'FontWeight','Bold','LabelFormat','none',...
    'LabelRotation','on','GLineWidth',.5,'GLineStyle','-','GAltitude',4,...
    'MLabelParallel',minlat-.01,'PLabelMeridian',minlon-.1);axis off; framem on; gridm off; mlabel on; plabel on;

surfm(rhoi.Lat_reg,rhoi.Lon_reg,(rhoi.di2d-rhoi.d02d)*1e9), shading interp
c=colorbar(ax2,'Position',[0.7 0.65 0.015 0.15]); ylabel(c, 'Residuals [E]'); 
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp, colormap(ax2,flipud(rwbcmap))
plotm(plat,plon,'k'),plotm(study_area(:,2),study_area(:,1),'k','LineWidth',2), 
plotm(transects2(1).lat,transects2(1).lon,'Color',[0.6 0.7 0.6],'LineWidth',2)
textm(transects2(1).lat(1)-0.1,transects2(1).lon(1)+0.1,'1','Color',[0 0 0],'FontSize',12,'FontWeight','bold')
set(gcf,'units','normalized','outerposition',[0.1 0.1 .5 .5])
print('-dpng','-r400','fig14')
print('-depsc','-painters','-r300','fig14')