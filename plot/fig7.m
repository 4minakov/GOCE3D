 %% Figure 7. Moho and Lithospheric cooling

%% Crustal thickness
Moho = -M_comb';
ff = scatteredInterpolant(Phi(:),Lam1(:),Moho(:));
fmoho = ff(lat2d,lon2d);
%
%ff = scatteredInterpolant(MLat(:),MLon1(:),rho_c1(:));
%rho_c_i = ff(lat2d,lon2d);
%rho_c_i = smoothn(rho_c_i,10);
%rho_c_i = rho_c_i*0+2850;
f1 = -fmoho.*(rho_c_i-rho_m);
f2 = -fmoho.^2.*(rho_c_i-rho_m);
f3 = -fmoho.^3.*(rho_c_i-rho_m);
% moho SH coefficients
cs1 = gsha(f1,method,grid,lmax);
cs2 = gsha(f2,method,grid,lmax);
cs3 = gsha(f3,method,grid,lmax);
sc1 = cs2sc(cs1);sc2 = cs2sc(cs2);sc3 = cs2sc(cs3);
% potential coefficients
K = 4*pi*R0^3./(M_E*(2*l+1));
Clm1 =        K.*sc1/R0^1;
Clm2 = Clm1 + K.*sc2/R0^2 .*(l+2) ;
Clm_moho_reg = Clm2 + K.*sc3/R0^3 .*(l+2).*(l+1);
%synthesis
[dg_moho, theRAD, lamRAD] = gshs_(Clm_moho_reg, 'quant', 'dg', ...
    'grid', 'block', 'gridsize', 2*lmax, 'height', 220e3, 'sub_wgs84', false);
%
[trr_moho, theRAD, lamRAD] = gshs_(Clm_moho_reg, 'quant', 'trr', ...
    'grid', 'block', 'gridsize', 2*lmax, 'height', 220e3, 'sub_wgs84', false);
%
dg_moho = dg_moho - mean(dg_moho(:));
trr_moho = trr_moho - mean(trr_moho(:));
%% thermal gravity correction
Clm_L90 = 0;
age_=age1';
ff=scatteredInterpolant(Phi(:),Lam1(:),age_(:));
age2=ff(Model.Lat,Model.Lon);
%
for i=1:Model.number_of_layers-1
     i
     Lupper = Model.(['l',num2str(i)]);
     Llower = Model.(['l',num2str(i+1)]);
     f = Lupper.bound-Llower.bound;
     depth = -mean(Lupper.bound(:)/2+Llower.bound(:)/2)
     %pause
     dens = -3400*1e-5*1400*...
         erf(depth./(2*sqrt(1e-6*age2*1e6*3600*24*365)) );
     dens = flipud(dens); 
     %dens = flipud(Lupper.dens); 
     %dens=dens-mean(dens(:));
     f1 = f.*dens;
     cs1 = gsha(f1,method,'block',lmax);
     cs1(1:4,1:4) = 0;
     sc1 = cs2sc(cs1);
     % potential coefficients
     K = 4*pi*R0^3./(M_E*(2*l+1));
     Clm1 = K.*sc1/R0^1;
     % correction for reference radius for each SH degree
     R2 = R0 + mean(mean((Lupper.bound)));
     Clm_L90 = Clm_L90 + Clm1.*(R2/R0).^l;     
end
 
 [dg_therm, theRAD, lamRAD] = gshs_(Clm_L90, 'quant', 'dg', ...
     'grid', 'block', 'gridsize', 2*lmax, 'height', 225e3, 'sub_wgs84', false);
 %
 [trr_therm, theRAD, lamRAD] = gshs_(Clm_L90, 'quant', 'trr', ...
     'grid', 'block', 'gridsize', 2*lmax, 'height', 225e3, 'sub_wgs84', false);
%%
figure('units','normalized','outerposition',[0 0 1 1])
colormap(roma)
%
subplot(2, 2, 1);
ax1=axesm('ortho','origin',[64 -22]); axis off
surfm(90-theRAD*180/pi, [lamRAD, lamRAD(1)]*180/pi, [dg_moho,dg_moho(:,1)]); title('a) \Delta g (Moho)'); 
plotm(study_area(:,2),study_area(:,1),'k','LineWidth',1.5),plotm([shorelines.Lat],[shorelines.Lon],'k')
c=colorbar(ax1,'FontSize',12,'Position',[0.45 0.6 0.01 0.26]);ylabel(c, '[mGal]');

subplot(2, 2, 2);
ax2=axesm('ortho','origin',[64 -22]); axis off
surfm(90-theRAD*180/pi, [lamRAD, lamRAD(1)]*180/pi, [trr_moho,trr_moho(:,1)]); title('b) T_{rr} (Moho)');
plotm(study_area(:,2),study_area(:,1),'k','LineWidth',1.5),plotm([shorelines.Lat],[shorelines.Lon],'k')
c=colorbar(ax2,'FontSize',12,'Position',[0.89 0.6 0.01 0.26]);ylabel(c, '[E]'); 

subplot(2, 2, 3);
ax3=axesm('ortho','origin',[64 -22]); axis off
surfm(90-theRAD*180/pi, [lamRAD, lamRAD(1)]*180/pi, [dg_therm,dg_therm(:,1)]); title('\Delta g (Thermal)'); 
plotm(study_area(:,2),study_area(:,1),'k','LineWidth',1.5),plotm([shorelines.Lat],[shorelines.Lon],'k')
c=colorbar(ax3,'FontSize',12,'Position',[0.45 0.15 0.01 0.26]);ylabel(c, '[mGal]');

subplot(2, 2, 4);
ax4=axesm('ortho','origin',[64 -22]); axis off
surfm(90-theRAD*180/pi, [lamRAD, lamRAD(1)]*180/pi, [trr_therm,trr_therm(:,1)]); title('T_{rr} (Thermal)'); 
plotm(study_area(:,2),study_area(:,1),'k','LineWidth',1.5),plotm([shorelines.Lat],[shorelines.Lon],'k')
c=colorbar(ax4,'FontSize',12,'Position',[0.89 0.15 0.01 0.26]);ylabel(c, '[E]');

set(gcf,'units','normalized','outerposition',[0 0 0.6 1])
print('-dpng','-r400','../fig/fig7')
print('-depsc','-r400','../fig/fig7')