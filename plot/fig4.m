%clear all
%close all

addpath(genpath('../data'))
addpath(genpath('../tools'))
load GSHHS_i
load XGM
load topography
load roma

% physical parameters
G = 6.67e-11;
minR=5971e3;
maxR=6371e3;
h_alt=220e3;
minlon=-59.5;
maxlon= 16.5;
minlat= 56.5;
maxlat= 82.5;
dlat=1;
dlon=4;
Lat_reg = maxlat:-dlat:minlat;
Lon_reg = minlon:dlon:maxlon;
nx=length(Lon_reg);
ny=length(Lat_reg);
ndata = nx*ny;
nr=10;
nmod = nx*ny*nr;
r = linspace(minR,maxR,nr);
dlon=dlon*pi/180;
dlat=dlat*pi/180;
dr=abs(r(2)-r(1));
[Lon3d,Lat3d,r3d]=meshgrid(Lon_reg,Lat_reg,r);
dV = dlon*dlat*dr;
[x3d,y3d,z3d]=sph2cart(Lon3d*pi/180,Lat3d*pi/180,r3d);
r_obs = maxR+h_alt;
[Xo,Yo,Zo]=sph2cart(Lon3d(:,:,1)*pi/180,Lat3d(:,:,1)*pi/180,r_obs);
%%
lam = Lon3d(:)*pi/180; lam_obs = Lon3d(:,:,1)*pi/180;
phi = Lat3d(:)*pi/180; phi_obs = Lat3d(:,:,1)*pi/180;
%
%datacovariancemodel_xyz %XYZ data covariance
sc = sqrt(cos(phi_obs(:))*cos(phi_obs(:))');
datacovariancemodel_rr %spherical data covariance
Cd = Cd_rr./sc;%single component data covariance
%
%% Green functions
tic
Gzz = zeros(ndata,nmod);Gxx=Gzz;Gyy=Gzz;Gxy=Gzz;Gzx=Gzz;Gzy=Gzz;ind=Gzz;
Grr=Gzz;L=Gzz;
%
dOm = r3d(:).^2 .* cos(pi/180*Lat3d(:))*dlat*dlon*dr;
%
for i=1:ndata
    %
    cospsi = sin(phi_obs(i))*sin(phi)+cos(phi_obs(i))*cos(phi).*cos(lam_obs(i)-lam);
    L(i,:)  = sqrt(r_obs.^2 + r3d(:).^2 - 2*r3d(:).*r_obs.*cospsi);
    %
    Grr(i,:) = r3d(:).^2.*cos(phi)./L(i,:)'.^5 .* (2*r_obs.^2 - r3d(:).^2 - 4*r_obs*r3d(:).*cospsi + 3*r3d(:).^2.*cospsi.^2);
    %
    rr = sqrt( (x3d(:)-Xo(i)).^2+(y3d(:)-Yo(i)).^2+(z3d(:)-Zo(i)).^2 );
    r2 = rr.*rr;
    r3 = r2.*rr;
    % Gravity gradients tensor Green's functions
    Gxx(i,:) = - 1./r3.*(1 - 3*(Xo(i)-x3d(:)).^2./r2);
    Gyy(i,:) = - 1./r3.*(1 - 3*(Yo(i)-y3d(:)).^2./r2);
    Gzz(i,:) = - 1./r3.*(1 - 3*(Zo(i)-z3d(:)).^2./r2);
    Gzy(i,:) =   1./r3.*(3*(Zo(i)-z3d(:)).*(Yo(i)-y3d(:))./r2);
    Gzx(i,:) =   1./r3.*(3*(Xo(i)-x3d(:)).*(Zo(i)-z3d(:))./r2);
    Gxy(i,:) =   1./r3.*(3*(Xo(i)-x3d(:)).*(Yo(i)-y3d(:))./r2);
end
toc
%%
%% model covariance
% synthetic data
%statistical model
%I  = 0.03; %radial correlation length
%Za = 0.95; 
%setup_name = ['synres','I',num2str(I),'Za',num2str(Za)];
nm = nmod;
lat1 = Lat3d*pi/180;
lon1 = Lon3d*pi/180;
r1  = r3d;
cm = 10^2*eye(nmod);
sigma_rho = 20^2;
lr = 200e3;%radial correlation length
w  = 2*pi/600e3;
Za = 0.94;%angular correlation parameter
tic
perc_complete = 0;
for k = 1:nmod
    cospsi = sin(lat1(k))*sin(lat1(:))+...
        cos(lat1(k))*cos(lat1(:)).*cos(lon1(k)-lon1(:));
    Cm_a = (1 - 2*Za*cospsi+Za^2).^(-1/2);Cm_a=Cm_a/max(Cm_a);
    % Radial correlation
    Cm_r = exp(-abs(r1(k)-r1(:))/lr).*cos(w*abs(r1(k)-r1(:)));
    %Cm_r = exp(-abs(r1(k)-r1(:))/lr);
    Cm_r=Cm_r/max(Cm_r);
    cm(k,:) = sigma_rho*Cm_a.*Cm_r;%local covariance matrix
    
    if round(k/nmod*100) > perc_complete
        disp([num2str(round(k/nmod*100)),'% complete..'])
        perc_complete = round(k/nmod*100);
    end
    
end
cm0 = cm;

if 1==1
    load drhoRND.mat
    % random realizations from covariance function
else
    if 1==1
        %Cholesky decomposition
        Nkl = nmod;
        LL = chol(cm0);
        LL = LL';
    else
        Nkl = 100;
        [V1,D1] = eigs(cm0,Nkl);
        D1=sqrt(D1);
        LL = V1*D1;
    end
    
    Q_ksi = randn(Nkl,1);
    ksi = LL*Q_ksi;
    m0  = 0*r3d;
    drho_rnd = m0+reshape(ksi,[ny nx nr]);
end
%%
figure, plot_spher2(drho_rnd,Lat_reg,Lon_reg,r)
colormap(roma),
ishore = find([shorelines.Lon] > minlon &...
    [shorelines.Lon] < maxlon & ...
    [shorelines.Lat] > minlat & ...
    [shorelines.Lat] < maxlat);
sh_lon = [shorelines.Lon];
sh_lat = [shorelines.Lat];
[sh_x,sh_y,sh_z]=sph2cart(sh_lon(ishore)*pi/180,sh_lat(ishore)*pi/180,sh_lat(ishore)*0+6375e3);
plot3(sh_x,sh_y,sh_z,'.k')
view(90,10), title('Input model')
%caxis([-50 50])
camlight('headlight'),%alpha(0.5)
c=colorbar('south');c.Position=[0.4 0.06 0.21 0.034];c.FontSize=12;
%%

%[Xm,Ym,Zm]=sph2cart( mean(Lon3d(:))*pi/180, mean(Lat3d(:))*pi/180, mean(r3d(:))    );
%drho_rnd =  -30 * exp(-((x3d-Xm).^2/300e3^2 + (y3d-Ym).^2/300e3^2 + (z3d-Zm).^2/100e3^2 ));
% forward problem
Trr = G*Grr*(drho_rnd(:).*dV);
Txx = G*Gxx*(drho_rnd(:).*dOm);
Tzz = G*Gzz*(drho_rnd(:).*dOm);
Tyy = G*Gyy*(drho_rnd(:).*dOm);
Txy = G*Gxy*(drho_rnd(:).*dOm);
Tzx = G*Gzx*(drho_rnd(:).*dOm);
Tzy = G*Gzy*(drho_rnd(:).*dOm);
%% LSQR model space
%model and data weighting
dS = (6371e3+220e3).^2 .* cos(pi/180*Lat3d(:,:,end))*dlat*dlon*dr;
Wm = (6371e3 - r3d + 220e3);%spherical weighting  - Trr
Wm = spdiags(Wm(:),0,nmod,nmod);
d_vec = [Trr(:)+0.1*std(Trr(:))*(1-2*rand(ndata,1))];%single component
ndata=length(Tzz(:));
Wd = 1./diag(Cd);
Wd = spdiags(Wd(:),0,ndata,ndata);%single component
%
lambda = 1e8;
rho0 = x3d*0;
d0 = d_vec*0;
GG = Grr;%single component Trr
A1 = G*GG*dV*Wm;%Gradient spher full
K = diag(diag(cm0,0));
AA = K*A1'*Wd;
A =   AA*A1 + lambda*speye(nmod,nmod);
d =  -AA*(d0-d_vec);
AA = [];
m_it = lsqr(A,d,1e-2,50);
m    = full(diag(Wm)).*m_it;
data_it = G*Grr*(m.*dV);
mLSQR=reshape(m,ny,nx,nr);
disp(['LSQR rms: data residual=',num2str(100*rms(d_vec-data_it)/rms(d_vec)),' % '...
    '  drho = ',num2str(rms(m)), ' kg m^-3' ])
%
figure,
plot_spher2(mLSQR,Lat_reg,Lon_reg,r), title('Model space (LSQR)'), colormap(roma), %caxis([-10 12])
ishore = find([shorelines.Lon] > minlon & [shorelines.Lon] < maxlon & ...
    [shorelines.Lat] > minlat & [shorelines.Lat] < maxlat);
sh_lon = [shorelines.Lon];
sh_lat = [shorelines.Lat];
[sh_x,sh_y,sh_z]=sph2cart(sh_lon(ishore)*pi/180,sh_lat(ishore)*pi/180,sh_lat(ishore)*0+6375e3);
plot3(sh_x,sh_y,sh_z,'.k')
view(90,10),camlight('headlight'),
c=colorbar('south');c.Position=[0.4 0.06 0.21 0.034];c.FontSize=12;
ylabel(c,'[kg m$^{-3}$]','Interpreter','latex','FontSize',14,'Position',[14 1])
%% inversion one-by-one data points
nd = ndata;%single component
d0 = d_vec+0.1*std(d_vec)*(1-2*rand(nd,1));
cd = diag(Cd);
wm=diag(Wm);
g  = G*GG*dV;
m0 = x3d(:)*0;
lambda=1e11;%spherical
cm = cm0;%spherical covariance
q = m0*0;
for k = 1:nd
    v   = d0(k) - (g(k,:).*wm')*(m0(:)./wm(:));%residual data value
    q   = cm(:,:)*(wm.*g(k,:)');
    a   = (wm'.*g(k,:))*q + lambda*cd(k);%updated hessian row
    m0  = m0 + wm(:).*(q*v/a);
    cm  = cm - q*q'/a;
    disp(['k=',num2str(k), ' data point, max=',num2str(max(abs(m0)))])
end
mPW = reshape(m0,ny,nx,nr);
data_it = G*Grr*(m0.*dV);
disp(['POINTWISE rms: data residual=',num2str(100*rms(d_vec-data_it)/rms(d_vec)),' % '...
    '  drho = ',num2str(rms(m0)), ' kg m^-3' ])

figure,
plot_spher2(mPW,Lat_reg,Lon_reg,r), title('Data space (pointwise)'), colormap(roma), %caxis([-10 12])
ishore = find([shorelines.Lon] > minlon & [shorelines.Lon] < maxlon & ...
    [shorelines.Lat] > minlat & [shorelines.Lat] < maxlat);
sh_lon = [shorelines.Lon];
sh_lat = [shorelines.Lat];
[sh_x,sh_y,sh_z]=sph2cart(sh_lon(ishore)*pi/180,sh_lat(ishore)*pi/180,sh_lat(ishore)*0+6375e3);
plot3(sh_x,sh_y,sh_z,'.k')
view(90,10),camlight('headlight'), %caxis([-10 10])
c=colorbar('south');c.Position=[0.4 0.06 0.21 0.034];c.FontSize=12;
ylabel(c,'[kg m$^{-3}$]','Interpreter','latex','FontSize',14,'Position',[14 1])
%% SVD full matrix inversion in data space
lambda = (1/length(diag(Wm)) * sum(1./diag(Wm).^2)).^-1;
A = G*Grr*Wm*dV;
hh = A*cm0*A' + lambda*Cd;
[U,S,V] = svd(hh);
s = diag(S);
%tol = max(size(hh)) * eps(norm(s,inf));
logtol = -21;
r1 = sum(log(s) > logtol);
V(:,r1:end) = [];
U(:,r1:end) = [];
s(r1:end) = [];
s = 1./s(:);
hh_inv = (V.*s.')*U';
cm1 = cm0 - cm0*A'*hh_inv*A*cm0;
m1 =  Wm*cm0*A'*hh_inv*d_vec;
mSVD = reshape(m1,ny,nx,nr);
data_it = G*Grr*(m1.*dV);
disp(['SVD rms: data residual=',num2str(100*rms(d_vec-data_it)/rms(d_vec)),' % '...
    '  drho = ',num2str(rms(m1)), ' kg m^-3' ])

figure,
plot_spher2(mSVD,Lat_reg,Lon_reg,r), title('Data space (SVD)'), colormap(roma), %caxis([-10 12])
ishore = find([shorelines.Lon] > minlon & [shorelines.Lon] < maxlon & ...
    [shorelines.Lat] > minlat & [shorelines.Lat] < maxlat);
sh_lon = [shorelines.Lon];
sh_lat = [shorelines.Lat];
[sh_x,sh_y,sh_z]=sph2cart(sh_lon(ishore)*pi/180,sh_lat(ishore)*pi/180,sh_lat(ishore)*0+6375e3);
plot3(sh_x,sh_y,sh_z,'.k')
view(90,10),camlight('headlight'), %caxis([-10 10])
c=colorbar('south');c.Position=[0.4 0.06 0.21 0.034];c.FontSize=12;
ylabel(c,'[kg m$^{-3}$]','Interpreter','latex','FontSize',14,'Position',[14 1])
%%

figure
subplot(221),plot_spher2(drho_rnd,Lat_reg,Lon_reg,r), title('a)    Input model     '), 
view(70,10),camlight('headlight'),axis tight
subplot(222),plot_spher2(mPW,Lat_reg,Lon_reg,r), title('b)    Data space (PW)     '), 
view(70,10),camlight('headlight'),
subplot(223),plot_spher2(mSVD,Lat_reg,Lon_reg,r), title('c)    Data space (SVD)    '), 
view(70,10),camlight('headlight'),
subplot(224),plot_spher2(mLSQR,Lat_reg,Lon_reg,r), title('d)    Model space (LSQR)'), 
view(70,10),camlight('headlight'),
colormap(roma)

set(gcf,'units','normalized','outerposition',[0 0 0.4 0.4])
print('-depsc','-r400','../fig/fig4')

if ~isfile('../data/DataSynth.mat')
    save ../data/DataSynth Lat_reg Lon_reg r drho_rnd mSVD mLSQR mPW Trr Tzz Txx Tzx Tzy Txy
end