addpath '..\tools'
addpath '..\calc'
addpath('C:\Program Files\gmt6\bin')


mod_real = load('real_data_TzzTxxTzx_2_results');

ndata = size(mod_real.Lat_reg,1).*size(mod_real.Lon_reg,1);
nmod = numel(mod_real.mPoint);
[ny,nx,nz] = size(mod_real.mPoint);
% G = 6.67e-11;
minR=5971e3;
maxR=6371e3;
h_alt=220e3; 
load Data_reg
Lat_reg = Lat_reg(1:1:end-2);
Lon_reg = Lon_reg(1:4:end-10);
Trr = 1e-9*Trr_reg(1:1:end-2,1:4:end-10);
dg_reg  = 1e-5*dg_reg(1:1:end-2,1:4:end-10);
nx=length(Lon_reg);
ny=length(Lat_reg);

dataECEF = load('GravTensDataRes_ECEF');
%titl = {'Vxx' 'Vxy' 'Vxz' 'Vyx' 'Vyy' 'Vyz' 'Vxz' 'Vyz' 'Vzz' };
[x2d0,y2d0]=pstereo(180/pi*dataECEF.theRAD2d,180/pi*dataECEF.lamRAD2d,70,0);
ff = scatteredInterpolant(x2d0(:),y2d0(:),dataECEF.T_res_xyz(:,1));
[lon2d,lat2d]=meshgrid(Lon_reg,Lat_reg);
[x2d,y2d]=pstereo(lat2d,lon2d,70,0);
Txx = ff(x2d,y2d); 
Txx = Txx-mean(Txx(:));

ff = scatteredInterpolant(x2d0(:),y2d0(:),dataECEF.T_res_xyz(:,2));
Txy = ff(x2d,y2d); 
Txy = Txy-mean(Txy(:));

ff = scatteredInterpolant(x2d0(:),y2d0(:),dataECEF.T_res_xyz(:,3));
Tzx = ff(x2d,y2d); 
Tzx = Tzx-mean(Tzx(:));

ff = scatteredInterpolant(x2d0(:),y2d0(:),dataECEF.T_res_xyz(:,5));
Tyy = ff(x2d,y2d); 
Tyy = Tyy-mean(Tyy(:));

ff = scatteredInterpolant(x2d0(:),y2d0(:),dataECEF.T_res_xyz(:,6));
Tzy = ff(x2d,y2d); 
Tzy = Tzy-mean(Tzy(:));

ff = scatteredInterpolant(x2d0(:),y2d0(:),dataECEF.T_res_xyz(:,9));
Tzz = ff(x2d,y2d); 
Tzz = Tzz-mean(Tzz(:));
%mod_real0 = mod_real; mod_real.mPoint=drho_rnd;

[Txx_rec,Tyy_rec,Tzz_rec,Tzy_rec,Tzx_rec,Txy_rec,Trr_rec] = forward_modeling(mod_real,mod_real.mPoint); 
%[Txx,Tyy,Tzz,Tzy,Tzx,Txy,Trr] = forward_modeling(mod_real0); 
%
% rmsVal = [
% rms(Trr-Trr_rec)/rms(Trr(:)),
% rms(Tzz-Tzz_rec)/rms(Tzz(:)),
% rms(Txx-Txx_rec)/rms(Txx(:)),
% rms(Tzx-Tzx_rec)/rms(Tzx(:))]'
%%
%addpath('C:\Program Files\gmt6\bin')
%load 'C:\Users\alexamin\Dropbox (UiO)\ESA2\Manuscript\REVISION2\fig\figData'
x = mod_real.Lon_reg;
y = mod_real.Lat_reg;
z = reshape(1e9*Trr_rec,ny,nx);
grdwrite2(x,y,z,'../data/Trr_rec_real.nc')
z = reshape(1e9*Tzz_rec,ny,nx);
grdwrite2(x,y,z,'../data/Tzz_rec_real.nc')
z = reshape(1e9*Txx_rec,ny,nx);
grdwrite2(x,y,z,'../data/Txx_rec_real.nc')
z = reshape(1e9*Tzx_rec,ny,nx);
grdwrite2(x,y,z,'../data/Tzx_rec_real.nc')
%
z = reshape(1e9*Trr,ny,nx);
grdwrite2(x,y,z,'../data/Trr_real.nc')
z = reshape(1e9*Tzz,ny,nx);
grdwrite2(x,y,z,'../data/Tzz_real.nc')
z = reshape(1e9*Txx,ny,nx);
grdwrite2(x,y,z,'../data/Txx_real.nc')
z = reshape(1e9*Tzx,ny,nx);
grdwrite2(x,y,z,'../data/Tzx_real.nc')
%
cd('../tools/'),pwd
system('fig18plot.bat')
cd('../plot/'),pwd
open('../fig/fig18gmt.pdf')