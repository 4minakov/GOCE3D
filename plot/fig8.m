%% Figure 8. Residual Anomaly and Seismic Tomography 
%% XGM
gfc_file = 'XGM2016.gfc';
[gsm, lmax, lmin, info] = parse_icgem(gfc_file, 'max_lm', 180); % read model
[field, lmax] = clm2sc(gsm, 'max_lm', lmax); % convert format of coeffs to /S|C\
field1=cs2sc(field); field1(1:3,1:3)=0;
Clm_xgm=field;
V_dg = gshs_(field1, 'quant', 'dg', 'grid', 'block', 'gridsize', 2*lmax, 'height', 225e3, 'sub_wgs84', 1);
%% GOCE gradients from grids
load GoceGrids% GOCE gradients from grids
NoL = 1801;
NoP = 900; 
LAT = reshape(T.data(:,2),NoL,NoP);
LON = reshape(T.data(:,1),NoL,NoP);
ZZA = reshape(T.data(:,6),NoL,NoP);
LAT1=LAT(1:end-1,:);
LON1=LON;LON1(LON<0)=LON1(LON<0)+360;LON1=LON1(1:end-1,:);
ZZA1=ZZA(1:end-1,:);
ff = scatteredInterpolant;
ff.Points=[LAT1(:),LON1(:)];
ff.Values=ZZA1(:);
[lamRAD2d,theRAD2d]=meshgrid(lamRAD,theRAD);
Trr_grd = ff(90-theRAD2d*180/pi,lamRAD2d*180/pi);
%%
dg_correction = dg_topo+dg_ice+dg_bath+dg_moho+dg_sed+dg_therm;
trr_correction = trr_topo+trr_ice+trr_bath+trr_moho+trr_sed+trr_therm;

dg_residual = V_dg - dg_correction;
trr_residual = Trr_grd - trr_correction;
%%

tomo = load('SL2013_dat');
%

clrmap = load('vik');
vik = interp1(1:size(clrmap.vik,1),clrmap.vik,linspace(1,size(clrmap.vik,1),12 ));

figure('units','normalized','outerposition',[0 0 1 1])
colormap(flipud(vik))
%
subplot(2, 2, 1);
ax1=axesm('ortho','origin',[64 -22]); axis off
surfm(90-theRAD*180/pi, [lamRAD, lamRAD(1)]*180/pi, [dg_residual,dg_residual(:,1)]); title('a) \Delta g (Residual)'); 
plotm(study_area(:,2),study_area(:,1),'k','LineWidth',1.5),plotm([shorelines.Lat],[shorelines.Lon],'k')
c=colorbar(ax1,'FontSize',12,'Position',[0.45 0.6 0.01 0.26]);ylabel(c, '[mGal]');

subplot(2, 2, 2);
ax2=axesm('ortho','origin',[64 -22]); axis off
surfm(90-theRAD*180/pi, [lamRAD, lamRAD(1)]*180/pi, [trr_residual,trr_residual(:,1)]); title('b) T_{rr} (Residual)');
plotm(study_area(:,2),study_area(:,1),'k','LineWidth',1.5),plotm([shorelines.Lat],[shorelines.Lon],'k')
c=colorbar(ax2,'FontSize',12,'Position',[0.89 0.6 0.01 0.26]);ylabel(c, '[E]'); 

subplot(2, 2, 3);
ax3=axesm('ortho','origin',[64 -22]); axis off
surfm(tomo.Lat, tomo.Lon, tomo.aaa(:,:,6)); title('\delta v_s (260 km)'); 
plotm(study_area(:,2),study_area(:,1),'k','LineWidth',1.5),plotm([shorelines.Lat],[shorelines.Lon],'k')
c=colorbar(ax3,'FontSize',12,'Position',[0.45 0.15 0.01 0.26]);ylabel(c, '[%]');

subplot(2, 2, 4);
ax4=axesm('ortho','origin',[64 -22]); axis off
surfm(tomo.Lat, tomo.Lon, tomo.aaa(:,:,2)); title('\delta v_s (80 km)'); 
plotm(study_area(:,2),study_area(:,1),'k','LineWidth',1.5),plotm([shorelines.Lat],[shorelines.Lon],'k')
c=colorbar(ax4,'FontSize',12,'Position',[0.89 0.15 0.01 0.26]);ylabel(c, ' [%]');

%set(gcf,'units','normalized','outerposition',[0 0 0.6 1])
print('-dpng','-r400','../fig/fig8')
print('-depsc','-r400','../fig/fig8')