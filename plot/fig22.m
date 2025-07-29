% figure 22
% Cross-plot drho-dVs
% calculated 2D histogram
pr1rhF = pr1rh+drhot1;
max_age = 80;
[age2d,depth2d]=meshgrid(transects2(1).age, tomoNA.depth);
age2d = interp1(transects2(1).dist,age2d',transects2(1).d1);age2d=age2d';
depth2d = interp1(transects2(1).dist,depth2d',transects2(1).d1);depth2d=depth2d';
xplot = dvs1([age2d<max_age & depth2d<250 & depth2d>50]);
yplot = pr1rhF([age2d<max_age & depth2d<250& depth2d>50])/3400*100;
hh = histogram2(xplot,yplot,20,'DisplayStyle','tile',...
    'Normalization','pdf','BinMethod','fd','EdgeAlpha',0.5, 'Visible','off');
xhh = hh.XBinEdges(1:end-1)+diff(hh.XBinEdges)/2;
yhh = hh.YBinEdges(1:end-1)+diff(hh.YBinEdges)/2;
%

figure, contourf(xhh,yhh, hh.Values'),shading interp
close(figure(1))
xlabel('\delta V_s [%]'), ylabel('\delta\rho [%]'), title('PDF')
colormap(flipud(bwcmap)),colorbar, hold on, plot(xhh, xhh*0.2,'r','LineWidth',2), axis([-10 5 -1.5 .7])
set(gcf,'units','normalized','outerposition',[0.1 0.1 .3 .5])
print('-dpng','-r400','../fig/fig22')
print('-depsc','-painters','-r300','../fig/fig22')