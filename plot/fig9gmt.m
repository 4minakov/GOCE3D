%% Figure 9. Enlarged figure for N Atlantic region. 
%
x = lamRAD*180/pi;
y = 90-theRAD*180/pi;
z = trr_bath+trr_topo+trr_ice;
grdwrite2(x,y,z,'../data/Trr_topofull.nc')
z = trr_sed;
grdwrite2(x,y,z,'../data/Trr_sed.nc')
z = trr_moho;
grdwrite2(x,y,z,'../data/Trr_moho.nc')
z = trr_therm;
grdwrite2(x,y,z,'../data/Trr_therm.nc')

cd('../tools/'),pwd
system('fig9plot.bat')
cd('../plot/'),pwd