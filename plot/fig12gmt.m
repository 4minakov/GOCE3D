% Plots result figures as in manuscript by
%
% Minakov, A. and Gaina, C.   
% Probabilistic linear inversion of satellite 
% gravity gradient data applied to the northeast Atlantic, 
% submitted to JGR Solid Earth, 2021JB021854RR
% 
% Last modified by alexamin@uio.no, 31/10/2021
%
% version v1.0

clear all
close all

addpath '..\tools'

addpath '..\data'

addpath('C:\Program Files\gmt6\bin')

load figData 
load rwbcmap
load bwcmap

addpath('C:\Program Files\gmt6\bin')

x = rhoi.Lon_reg;
y = rhoi.Lat_reg;
z = 1e9*sqrt(rhoi.cd2d);
grdwrite2(x,y,z,'../data/trr_cd.nc')
%%
system('del gmt.history')

G = gmt('read -Tg ../data/trr_cd.nc'); P=[];
%gmt('makecpt -Cgray12 -T0.4/1.1/.1 > graycolors');
%gmt('makecpt -Cgray12 -T0.4/1.1/.1 > graycolors');
gmt('grd2cpt -Cgray > graycolors',G)
%
%gmt('grdcontour -R-55/20/55/82 -Js-20/90/12c/60 -B30/10 -P -E150 -Ngraycolors -Cgraycolors -K > fig12new.eps', G);
gmt('grdimage -R-55/20/55/82 -Js-20/90/12c/60 -B30/10 -P -E150 -Cgraycolors -K > ../fig/fig12new.eps', G);
gmt('psscale -R -J -Cgraycolors -O -K -B0.1x -By+lE >> ../fig/fig12new.eps');
gmt('pscoast','-Di -W0.5p -A1000 -N1/0.5p -J -R -W0.5p,yellow -O >> ../fig/fig12new.eps')
gmt('psconvert','../fig/fig12new.eps -Tf -P -A ')
open('../fig/fig12new.pdf')