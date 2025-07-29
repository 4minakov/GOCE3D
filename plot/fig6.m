%% Figure 6. Topography and Sediments 
%% bathymetry
Lam1=Lam;Lam1(Lam<0)=Lam1(Lam<0)+360;
topo1=topoi';
ff = scatteredInterpolant(Phi(:),Lam1(:),topo1(:));
f = ff(lat2d,lon2d);
f(f>0)=0;
%rho_c = f*0+2850;
MLon1=MLon;MLon1(MLon<0)=MLon1(MLon<0)+360;
rho_c1 = rho_c';
ff = scatteredInterpolant(MLat(:),MLon1(:),rho_c1(:));
rho_c_i = ff(lat2d,lon2d);
rho_c_i = smoothn(rho_c_i,10);
%
f1 = -f.*(rho_w-rho_c_i);
f2 = -f.^2.*(rho_w-rho_c_i);
f3 = -f.^3.*(rho_w-rho_c_i);
% bathymetry SH coefficients
cs1 = gsha(f1,method,grid,lmax);
cs2 = gsha(f2,method,grid,lmax);
cs3 = gsha(f3,method,grid,lmax);
%cs(1:4,1:4) = 0;
sc1 = cs2sc(cs1);sc2 = cs2sc(cs2);sc3 = cs2sc(cs3);
% potential coefficients
K = 4*pi*R0^3./(M_E*(2*l+1));
Clm1 =        K.*sc1/R0^1;
Clm2 = Clm1 + K.*sc2/R0^2 .*(l+2) ;
Clm_bath = Clm2 + K.*sc3/R0^3 .*(l+2).*(l+1);
%synthesis
[dg_bath, theRAD, lamRAD] = gshs_(Clm_bath, 'quant', 'dg', ...
    'grid', 'block', 'gridsize', 2*lmax, 'height', 220e3, 'sub_wgs84', false);
%
[trr_bath, theRAD, lamRAD] = gshs_(Clm_bath, 'quant', 'trr', ...
    'grid', 'block', 'gridsize', 2*lmax, 'height', 220e3, 'sub_wgs84', false);
%
dg_bath = dg_bath-mean(dg_bath(:));
trr_bath = trr_bath-mean(trr_bath(:));
%% topography
Lam1=Lam;Lam1(Lam<0)=Lam1(Lam<0)+360;
topo1=topoi';
ff = scatteredInterpolant(Phi(:),Lam1(:),topo1(:));
f = ff(lat2d,lon2d);
f(f<0)=0;
f1 = f.*rho_m;
f2 = f.^2.*rho_m;
f3 = f.^3.*rho_m;
% topography SH coefficients
cs1 = gsha(f1,method,grid,lmax);
cs2 = gsha(f2,method,grid,lmax);
cs3 = gsha(f3,method,grid,lmax);
sc1 = cs2sc(cs1);sc2 = cs2sc(cs2);sc3 = cs2sc(cs3);
% potential coefficients
K = 4*pi*R0^3./(M_E*(2*l+1));
Clm1 =        K.*sc1/R0^1;
Clm2 = Clm1 + K.*sc2/R0^2 .*(l+2) ;
Clm_topo = Clm2 + K.*sc3/R0^3 .*(l+2).*(l+1);
%synthesis
[dg_topo, theRAD, lamRAD] = gshs_(Clm_topo, 'quant', 'dg', ...
    'grid', 'block', 'gridsize', 2*lmax, 'height', 220e3, 'sub_wgs84', false);
%
[trr_topo, theRAD, lamRAD] = gshs_(Clm_topo, 'quant', 'trr', ...
    'grid', 'block', 'gridsize', 2*lmax, 'height', 220e3, 'sub_wgs84', false);
%% ice
rho_ice = 970;
ice_top = topoi'; ice_top(Ice<=0)=0;
ice_bot = ice_top-Ice';
ff = scatteredInterpolant(Phi(:),Lam1(:),ice_top(:));
f_t = ff(lat2d,lon2d);
ff = scatteredInterpolant(Phi(:),Lam1(:),ice_bot(:));
f_b = ff(lat2d,lon2d);
f1 = (f_t-f_b).*(rho_ice-rho_m);
f2 = (f_t.^2-f_b.^2).*(rho_ice-rho_m);
f3 = (f_t.^3-f_b.^3).*(rho_ice-rho_m);
% topography SH coefficients
cs1 = gsha(f1,method,grid,lmax);
cs2 = gsha(f2,method,grid,lmax);
cs3 = gsha(f3,method,grid,lmax);
sc1 = cs2sc(cs1);sc2 = cs2sc(cs2);sc3 = cs2sc(cs3);
% potential coefficients
K = 4*pi*R0^3./(M_E*(2*l+1));
Clm1 =        K.*sc1/R0^1;
Clm2 = Clm1 + K.*sc2/R0^2 .*(l+2) ;
Clm_ice = Clm2 + K.*sc3/R0^3 .*(l+2).*(l+1);
%synthesis
[dg_ice, theRAD, lamRAD] = gshs_(Clm_ice, 'quant', 'dg', ...
    'grid', 'block', 'gridsize', 2*lmax, 'height', 220e3, 'sub_wgs84', false);
[trr_ice, theRAD, lamRAD] = gshs_(Clm_ice, 'quant', 'trr', ...
    'grid', 'block', 'gridsize', 2*lmax, 'height', 220e3, 'sub_wgs84', false);
%% sediments with depth-dependent density
dens_prof = 'constant';%linear, quadratic
%dens_prof = 'linear';
%dens_prof = 'quadratic';
outfile = ['../data/trr_sed_',dens_prof];
if strcmp(dens_prof,'constant')
    a = 0; b = 0; rho_s0 = rho_s-rho_c0;
elseif strcmp(dens_prof,'linear')
    a = -6e-5; b = 0; rho_s0 = 2.2e3-rho_c0;
elseif strcmp(dens_prof,'quadratic')
    a = -9e-5; b = 2.5e-9; rho_s0 = 2.2e3-rho_c0;
end
%
sed_top = topoi'; 
sed_bot = double(topoi'-sed_comb');
ff = scatteredInterpolant(Phi(:),Lam1(:),sed_top(:));
f_t = ff(lat2d,lon2d);
ff = scatteredInterpolant(Phi(:),Lam1(:),sed_bot(:));
f_b = ff(lat2d,lon2d);
f1 = (f_t-f_b).*rho_s0;
f2 = (f_t.^2-f_b.^2).*rho_s0;
f3 = (f_t.^3-f_b.^3).*rho_s0;
% topography SH coefficients
cs1 = gsha(f1,method,grid,lmax);
cs2 = gsha(f2,method,grid,lmax);
cs3 = gsha(f3,method,grid,lmax);
sc1 = cs2sc(cs1);sc2 = cs2sc(cs2);sc3 = cs2sc(cs3);
% potential coefficients
K = 4*pi*R0^3./(M_E*(2*l+1));
Clm1 =        K.*sc1/R0^1;
Clm2 = Clm1 + K.*sc2/2/R0^2 .*((l+2)-a*R0) ;
Clm_sed = Clm2 + K.*sc3/6/R0^3 .*((l+2).*(l+1) - 2*(l+2)*a*R0 + 2*b*R0^2);
%synthesis
[dg_sed, theRAD, lamRAD] = gshs_(Clm_sed, 'quant', 'dg', ...
    'grid', 'block', 'gridsize', 2*lmax, 'height', 220e3, 'sub_wgs84', false);
[trr_sed, theRAD, lamRAD] = gshs_(Clm_sed, 'quant', 'trr', ...
    'grid', 'block', 'gridsize', 2*lmax, 'height', 220e3, 'sub_wgs84', false);
%%

clrmap = load('roma');
roma = interp1(1:size(clrmap.roma,1),clrmap.roma,linspace(1,size(clrmap.roma,1),12 ));

dg_topo_total = dg_topo+dg_ice+dg_bath;
trr_topo_total = trr_topo+trr_ice+trr_bath;

figure('units','normalized','outerposition',[0 0 1 1])
colormap(roma)
%
subplot(2, 2, 1);
ax1=axesm('ortho','origin',[64 -22]); axis off
surfm(90-theRAD*180/pi, [lamRAD, lamRAD(1)]*180/pi, [dg_topo_total,dg_topo_total(:,1)]); title('a) \Delta g (Topography)'); 
plotm(study_area(:,2),study_area(:,1),'k','LineWidth',1.5),plotm([shorelines.Lat],[shorelines.Lon],'k')
c=colorbar(ax1,'FontSize',12,'Position',[0.45 0.6 0.01 0.26]);ylabel(c, '[mGal]');

subplot(2, 2, 2);
ax2=axesm('ortho','origin',[64 -22]); axis off
surfm(90-theRAD*180/pi, [lamRAD, lamRAD(1)]*180/pi, [trr_topo_total,trr_topo_total(:,1)]); title('b) T_{rr} (topography)');
plotm(study_area(:,2),study_area(:,1),'k','LineWidth',1.5),plotm([shorelines.Lat],[shorelines.Lon],'k')
c=colorbar(ax2,'FontSize',12,'Position',[0.89 0.6 0.01 0.26]);ylabel(c, '[E]'); 

subplot(2, 2, 3);
ax3=axesm('ortho','origin',[64 -22]); axis off
surfm(90-theRAD*180/pi, [lamRAD, lamRAD(1)]*180/pi, [dg_sed,dg_sed(:,1)]); title('\Delta g (sediments)'); 
plotm(study_area(:,2),study_area(:,1),'k','LineWidth',1.5),plotm([shorelines.Lat],[shorelines.Lon],'k')
c=colorbar(ax3,'FontSize',12,'Position',[0.45 0.15 0.01 0.26]);ylabel(c, '[mGal]');

subplot(2, 2, 4);
ax4=axesm('ortho','origin',[64 -22]); axis off
surfm(90-theRAD*180/pi, [lamRAD, lamRAD(1)]*180/pi, [trr_sed,trr_sed(:,1)]); title('T_{rr} (sediments)'); 
plotm(study_area(:,2),study_area(:,1),'k','LineWidth',1.5),plotm([shorelines.Lat],[shorelines.Lon],'k')
c=colorbar(ax4,'FontSize',12,'Position',[0.89 0.15 0.01 0.26]);ylabel(c, '[E]');

set(gcf,'units','normalized','outerposition',[0 0 0.6 1])
print('-dpng','-r400','../fig/fig6')
print('-depsc','-r400','../fig/fig6')