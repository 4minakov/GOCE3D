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

% clear all
% close all

addpath '../tools'
addpath '../data'
addpath('C:\Program Files\gmt6\bin')

load figData 
load rwbcmap
load bwcmap

aa2 = load('trr_sed_quadratic');
aa0 = load('trr_sed_constant');
bb  = load('ForwModGrd');

x = bb.lamRAD*180/pi;
y = 90-bb.theRAD*180/pi;
z = (aa0.trr_sed-aa2.trr_sed);

grdwrite2(x,y,z,'../data/trr_sed02_res.nc')
%%
system('del gmt.history')

G = gmt('read -Tg ../data/trr_sed02_res.nc'); P=[];

%gmt('grd2cpt -Cgray > graycolors',G)
%gmt('grd2cpt -Croma > romacolors',G)
gmt('makecpt -Croma -Ic -T-0.5/0.5/.1 > colors');
%
%gmt('grdcontour -R-55/20/55/82 -Js-20/90/12c/60 -B30/10 -P -E150 -Ngraycolors -Cgraycolors -K > fig12new.eps', G);
gmt('grdimage -R-55/20/55/82 -Js-20/90/12c/60 -B30/10 -P -E150 -Ccolors -K > ../fig/fig10gmt.eps', G);
gmt('psscale -R -J -Ccolors -O -K -B0.1x -By+lE >> ../fig/fig10gmt.eps');
gmt('pscoast','-Di -W0.5p -A1000 -N1/0.5p -J -R -O >> ../fig/fig10gmt.eps')
gmt('psconvert','../fig/fig10gmt.eps -Tf -P -A ')
open('../fig/fig10gmt.pdf')