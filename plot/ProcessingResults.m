lmax = 180;

minlat = 55; maxlat = 84; nlat=maxlat-minlat;
minlon = -60; maxlon = 30; nlon=maxlon-minlon;

study_area = [linspace(minlon,maxlon,100)',maxlat*ones(100,1);...
    linspace(maxlon,minlon,100)',minlat*ones(100,1);...
    minlon,maxlat];
%
ny0=1081; nx0=2161;
lontop = linspace(-180,180,nx0); lattop = linspace(  90,-90,ny0);
[Phi,Lam] = meshgrid(lattop,lontop);
%
d = age; 
d(age<=75&age>0)       = -2652 - 324*sqrt(age(age<=75&age>0));
d(age>75&age<=160) = -5028 - 5.26*d(age>75&age<=160) + 250*sin((age(age>75&age<=160) - 75)/30);
d(age>160)=-5750;
d(age<0)=NaN;
w = d + 5750; 
w(isnan(w))=0; w(w<0)=0;
%
rho_c=Vp;
w1 = interp2(Phi,Lam,w',MLat,MLon); 
rho_c(w1>0) = 350 + 385*Vp(w1>0);
rho_c(w1==0) = 590 + 346*Vp(w1==0);
rho_c = rho_c';
%
rho_c0 = 2850;
rho_m = 3400;
rho_w = 1030;
rho_s = 2400;
R0 = 6371e3;%Earth mean radius
rho_E = 5510; %mean Earth density
M_E = 4/3*pi*R0.^3*rho_E;
method = 'aq';
grid = 'block';
l = (0:lmax)';
l = repmat(l,1,2*lmax+1);
% coordinate mesh
lon = linspace(0.5, 359.5, 720);
lat = linspace(89.5, -89.5, 360); 
[lon2d,lat2d]=meshgrid(lon,lat);

ii = find(tomoNA.depth==160);%depth slice
dvs_int = tomoNA.vsh(:,:,ii)';

ii_study =  find(Phi>=minlat & Phi<=maxlat & Lam>=minlon & Lam<=maxlon ) ;
Phi1=Phi(lontop>=minlon & lontop<=maxlon,lattop>=minlat & lattop<=maxlat);
Lam1=Lam(lontop>=minlon & lontop<=maxlon,lattop>=minlat & lattop<=maxlat);
topo1=topoi';
topo1=topo1(lontop>=minlon & lontop<=maxlon,lattop>=minlat & lattop<=maxlat);

jj=10;

drhoz = rhoi.m3d(:,:,jj);
[rhoi.Lat_reg2d,rhoi.Lon_reg2d]=ndgrid(rhoi.Lat_reg(:),rhoi.Lon_reg(:));
FF = scatteredInterpolant(rhoi.Lat_reg2d(:),rhoi.Lon_reg2d(:), drhoz(:));
FF.ExtrapolationMethod='none';
rho_int = FF(Phi1, Lam1);


