%% Figure 11. Residual Trr anomaly.
addpath('C:\Program Files\gmt6\bin')
z = trr_residual;
grdwrite2(x,y,z,'../data/Trr_residual.nc')

G = gmt('read -Tg ../data/Trr_residual.nc'); P=[];
gmt('makecpt -Cvik -Ic -T-3/3/.5 -Z > colors');
%
gmt('grdimage -R-55/20/55/82 -Js-20/90/12c/60 -B30/10 -P -E150 -Ccolors -K > ../fig/fig11gmt.eps', G);
gmt('psscale -R -J -Ccolors -O -K -Bx -By+lE >> ../fig/fig11gmt.eps');
gmt('pscoast','-Di -W0.5p -A1000 -N1/0.5p -J -R -O >> ../fig/fig11gmt.eps')
gmt('psconvert','../fig/fig11gmt.eps -Tf -P -A ')
open('../fig/fig11gmt.pdf')