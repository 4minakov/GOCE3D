% Plots result figures as in manuscript by
%
% Minakov, A. and Gaina, C.   
% Probabilistic linear inversion of satellite 
% gravity gradient data applied to the northeast Atlantic, 
% submitted to JGR Solid Earth, 2021JB021854RR
% 
% Last modified by alexamin@uio.no, 31/10/2021
%
% version v1.0m

addpath '..\tools'

addpath '..\data'

addpath('C:\Program Files\gmt6\bin')

load figData 
load rwbcmap
load bwcmap

addpath('C:\Program Files\gmt6\bin')

[Txx,Tyy,Tzz,Tzy,Tzx,Txy,Trr] = forward_modeling(rhoi,rhoi.m3d); 
%%
x = rhoi.Lon_reg;
y = rhoi.Lat_reg;
z = 1e9*(reshape(Trr,size(rhoi.d02d))-rhoi.d02d);
grdwrite2(x,y,z,'../data/trr_inv_res.nc')
%%
system('del gmt.history')

G = gmt('read -Tg ../data/trr_inv_res.nc'); P=[];

gmt('makecpt -Cvik -Ic -T-1/1/.2 -Z > colors');
%
%gmt('grdcontour -R-55/20/55/82 -Js-20/90/12c/60 -B30/10 -P -E150 -Ngraycolors -Cgraycolors -K > fig12new.eps', G);
gmt('grdimage -R-55/20/55/82 -Js-20/90/12c/60 -B30/10 -P -E150 -Ccolors -K > ../fig/fig14gmt.eps', G);
gmt('psscale -R -J -Ccolors -O -K -B0.5x -By+lE >> ../fig/fig14gmt.eps');
gmt('pscoast','-Di -W0.5p -A1000 -N1/0.5p -J -R -O >> ../fig/fig14gmt.eps')
gmt('psconvert','../fig/fig14gmt.eps -Tf -P -A ')
open('../fig/fig14gmt.pdf')