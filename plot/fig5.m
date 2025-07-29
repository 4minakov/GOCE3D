%clear all
%close all

addpath(genpath('../data'))
addpath(genpath('../tools'))
addpath('C:\Program Files\gmt6\bin')

load GSHHS_i
load XGM
load topography
load roma
mod_syn = load('DataSynth');

ndata = size(mod_syn.Lat_reg,1).*size(mod_syn.Lon_reg,1);
nmod = numel(mod_syn.mPW);
[ny,nx,nz] = size(mod_syn.mPW);
% G = 6.67e-11;
minR=5971e3;
maxR=6371e3;
h_alt=220e3; 

[Txx_rec,Tyy_rec,Tzz_rec,Tzy_rec,Tzx_rec,Txy_rec,Trr_rec] = forward_modeling(mod_syn,mod_syn.mPW); 
[Txx,Tyy,Tzz,Tzy,Tzx,Txy,Trr] = forward_modeling(mod_syn,mod_syn.drho_rnd); 

%%

x = mod_syn.Lon_reg;
y = mod_syn.Lat_reg;
z = reshape(1e9*Trr_rec,ny,nx);
grdwrite2(x,y,z,'../data/Trr_rec_syn.nc')
z = reshape(1e9*Tzz_rec,ny,nx);
grdwrite2(x,y,z,'../data/Tzz_rec_syn.nc')
z = reshape(1e9*Txx_rec,ny,nx);
grdwrite2(x,y,z,'../data/Txx_rec_syn.nc')
z = reshape(1e9*Tzx_rec,ny,nx);
grdwrite2(x,y,z,'../data/Tzx_rec_syn.nc')
%
z = reshape(1e9*Trr,ny,nx);
grdwrite2(x,y,z,'../data/Trr_syn.nc')
z = reshape(1e9*Tzz,ny,nx);
grdwrite2(x,y,z,'../data/Tzz_syn.nc')
z = reshape(1e9*Txx,ny,nx);
grdwrite2(x,y,z,'../data/Txx_syn.nc')
z = reshape(1e9*Tzx,ny,nx);
grdwrite2(x,y,z,'../data/Tzx_syn.nc')
%
cd('../tools/'),pwd
system('fig5plot.bat')
cd('..'),pwd

%
disp(['RMS misfit (%) : ', ...
' Trr = ', num2str(ceil(100*rms(Trr-Trr_rec)/rms(Trr(:)))),...
' Tzz = ', num2str(ceil(100*rms(Tzz-Tzz_rec)/rms(Tzz(:)))),...
' Txx = ', num2str(ceil(100*rms(Txx-Txx_rec)/rms(Txx(:)))),...
' Tzx = ', num2str(ceil(100*rms(Tzx-Tzx_rec)/rms(Tzx(:))))  ])
%
% gmt('makecpt -Croma -Ic -T-1/1/.2 -Z > colors');
% %
% G=gmt('read -Tg Trr_syn.nc');P=[];
% P=gmt('grdimage -R-56/14/57/80 -Js-20/90/12c/60 -B30/10 -P -E150 -Ccolors -K', G);
% P=gmt('psscale -R -J -Ccolors -O');
% Irr=gmt('psconvert','-TG -P -E600 -A', P);
% %
% G=gmt('read -Tg Tzz_syn.nc');P=[];
% P=gmt('grdimage -R-56/14/57/80 -Js-20/90/12c/60 -B30/10 -P -E150 -Ccolors -K', G);
% P=gmt('psscale -R -J -Ccolors -O');
% Izz=gmt('psconvert','-TG -P -E600 -A', P);
% %
% G=gmt('read -Tg Txx_syn.nc');P=[];
% P=gmt('grdimage -R-56/14/57/80 -Js-20/90/12c/60 -B30/10 -P -E150 -Ccolors -K', G);
% P=gmt('psscale -R -J -Ccolors -O');
% Ixx=gmt('psconvert','-TG -P -E600 -A', P);
% %
% G=gmt('read -Tg Tzx_syn.nc');P=[];
% P=gmt('grdimage -R-56/14/57/80 -Js-20/90/12c/60 -B30/10 -P -E150 -Ccolors -K', G);
% P=gmt('psscale -R -J -Ccolors -O');
% Izx=gmt('psconvert','-TG -P -E600 -A', P);
% %
% G=gmt('read -Tg Trr_rec_syn.nc');P=[];
% P=gmt('grdimage -R-56/14/57/80 -Js-20/90/12c/60 -B30/10 -P -E150 -Ccolors -K', G);
% P=gmt('psscale -R -J -Ccolors -O');
% IrrR=gmt('psconvert','-TG -P -E600 -A', P);
% %
% G=gmt('read -Tg Tzz_rec_syn.nc');P=[];
% P=gmt('grdimage -R-56/14/57/80 -Js-20/90/12c/60 -B30/10 -P -E150 -Ccolors -K', G);
% P=gmt('psscale -R -J -Ccolors -O');
% IzzR=gmt('psconvert','-TG -P -E600 -A', P);
% %
% G=gmt('read -Tg Txx_rec_syn.nc');P=[];
% P=gmt('grdimage -R-56/14/57/80 -Js-20/90/12c/60 -B30/10 -P -E150 -Ccolors -K', G);
% P=gmt('psscale -R -J -Ccolors -O');
% IxxR=gmt('psconvert','-TG -P -E600 -A', P);
% %
% G=gmt('read -Tg Tzx_rec_syn.nc');P=[];
% P=gmt('grdimage -R-56/14/57/80 -Js-20/90/12c/60 -B30/10 -P -E150 -Ccolors -K', G);
% P=gmt('psscale -R -J -Ccolors -O');
% IzxR=gmt('psconvert','-TG -P -E600 -A', P);
% %
% figure(1); clf
% subplot(421)
% h=imshow(Irr.image); 
% set (h, 'AlphaData', Irr.alpha)
% tt=title('(a)           T_{rr}');tt.Units='normalized';
% tt.Position(1)=0; tt.HorizontalAlignment='left';
% subplot(422),
% h=imshow(IrrR.image); 
% set (h, 'AlphaData', IrrR.alpha)
% tt=title('(b)');tt.Units='normalized';
% tt.Position(1)=0; tt.HorizontalAlignment='left';
% %
% subplot(423)
% h=imshow(Izz.image); 
% set (h, 'AlphaData', Izz.alpha)
% tt=title('(c)           T_{zz}');tt.Units='normalized';
% tt.Position(1)=0; tt.HorizontalAlignment='left';
% subplot(424),
% h=imshow(IzzR.image); 
% set (h, 'AlphaData', IzzR.alpha)
% tt=title('(d) ');tt.Units='normalized';
% tt.Position(1)=0; tt.HorizontalAlignment='left';
% %
% subplot(425)
% h=imshow(Ixx.image); 
% set (h, 'AlphaData', Ixx.alpha)
% tt=title('(e)           T_{xx}');tt.Units='normalized';
% tt.Position(1)=0; tt.HorizontalAlignment='left';
% subplot(426),
% h=imshow(IxxR.image); 
% set (h, 'AlphaData', IxxR.alpha)
% tt=title('(f) ');tt.Units='normalized';
% tt.Position(1)=0; tt.HorizontalAlignment='left';
% %
% subplot(427)
% h=imshow(Izx.image); 
% set (h, 'AlphaData', Izx.alpha)
% tt=title('(g)           T_{zx}');tt.Units='normalized';
% tt.Position(1)=0; tt.HorizontalAlignment='left';
% subplot(428),
% h=imshow(IzxR.image); 
% set (h, 'AlphaData', IzxR.alpha)
% tt=title('(h) ');tt.Units='normalized';
% tt.Position(1)=0; tt.HorizontalAlignment='left';
% 
% set(gcf,'Units','normalized','OuterPosition',[0.1 0.1 .25 0.75])
% print('-depsc','-painters','-r400','fig5')


