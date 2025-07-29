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

addpath 'C:\Users\alexamin\Dropbox (UiO)\ESA2\Manuscript\REVISION1\tools'
addpath 'C:\Users\alexamin\Dropbox (UiO)\ESA2\Manuscript\REVISION2\calc'
addpath 'C:\Users\alexamin\Dropbox (UiO)\ESA2\Manuscript\REVISION2\data'
addpath 'C:\Users\alexamin\Dropbox (UiO)\ESA2\Manuscript\REVISION2\fig'
addpath('C:\Program Files\gmt6\bin')

load figData 
load rwbcmap
load bwcmap

addpath('C:\Program Files\gmt6\bin')
load 'C:\Users\alexamin\Dropbox (UiO)\ESA2\Manuscript\REVISION2\fig\figData'

aa2 = load('trr_sed_quadratic');
aa0 = load('trr_sed_constant');
bb  = load('ForwModGrd');

x = rhoi.Lon_reg;
y = rhoi.Lat_reg;
z = 1e9*(rhoi.di2d-rhoi.d02d);
grdwrite2(x,y,z,'trr_inv_res.nc')
%%
system('del gmt.history')

G = gmt('read -Tg trr_inv_res.nc'); P=[];

%gmt('grd2cpt -Cgray > graycolors',G)
%gmt('grd2cpt -Croma > romacolors',G)
gmt('makecpt -Cvik -Ic -T-1/1/.2 -Z > romacolors');
%
%gmt('grdcontour -R-55/20/55/82 -Js-20/90/12c/60 -B30/10 -P -E150 -Ngraycolors -Cgraycolors -K > fig12new.eps', G);
gmt('grdimage -R-55/20/55/82 -Js-20/90/12c/60 -B30/10 -P -E150 -Cromacolors -K > fig14new.eps', G);
gmt('psscale -R -J -Cromacolors -O -K -B0.5x -By+lE >> fig14new.eps');
gmt('pscoast','-Di -W0.5p -A1000 -N1/0.5p -J -R -O >> fig14new.eps')
gmt('psconvert','fig14new.eps -Tf -P -A ')
open('fig14new.pdf')