% figure 15
% Density variation at 160 km  
figure,
ax2=axesm('lambertstd','MapLatLimit',[minlat maxlat-2],'MapLonLimit',[minlon maxlon-8],...
    'MLineLocation', 10, 'PLineLocation', 5,...
    'FontSize', 12, 'FontWeight','Bold','LabelFormat','none',...
    'LabelRotation','on','GLineWidth',.5,'GLineStyle','-','GAltitude',4,...
    'MLabelParallel',minlat-.01,'PLabelMeridian',minlon-.1);axis off; framem on; gridm off; mlabel on; plabel on;
surfm(Phi1,Lam1,rho_int,topo1*1e-5-100);
camlight, material dull, lighting Gouraud, shading interp, alpha(1.)
c=colorbar(ax2,'Position',[0.7 0.65 0.015 0.15]); ylabel(c, 'Density [kg m^{-3}]'); 
hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp, colormap(ax2,flipud(rwbcmap))
plotm(plat,plon,'k'),plotm(study_area(:,2),study_area(:,1),'k','LineWidth',2), %caxis([-20 20])
plotm(transects2(1).lat,transects2(1).lon,'Color',[0.6 0.7 0.6],'LineWidth',2)
textm(transects2(1).lat(1)-0.1,transects2(1).lon(1)+0.1,'1','Color',[0 0 0],'FontSize',12,'FontWeight','bold')
for i=2:2:28
    plotm(cot.(['l',num2str(i)])(:,2),cot.(['l',num2str(i)])(:,1), 'Color',[.3 .3 .3],'LineWidth',1)
end
plotm([volc_centr.Lat],[volc_centr.Lon],'y')
title(['Depth ',num2str(fix(6371-rhoi.r(jj)/1000)),' km'])
set(gcf,'units','normalized','outerposition',[0.1 0.1 .5 .5])
print('-dpng','-r400','../fig/fig15')
print('-depsc','-painters','-r300','../fig/fig15')
