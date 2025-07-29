%% data covariance matrix
phi0 =  90;%true latitude
lam0 = -0;%zero longitude
mkrig_var = load('test1_var.out');
[Xkrig,Ykrig]= pstereo(mkrig_var(:,2),mkrig_var(:,1),phi0,lam0);
FKRIG = scatteredInterpolant(Xkrig,Ykrig,mkrig_var(:,3));
FKRIG.Method = 'natural';
FKRIG.ExtrapolationMethod='none';
[Xgrd,Ygrd]=pstereo(phi_obs*180/pi,lam_obs*180/pi,phi0,lam0);
mvar_grd = FKRIG(Xgrd,Ygrd);
mvar_grd(isnan(mvar_grd))=0;
Cms = 1e6*diag(mvar_grd(:));
%mstd     = sqrt(mvar_grd);
drho_s   = 450;
mvar_grd =FKRIG(Xgrd,Ygrd);
%figure,
%contourf(lam_obs*180/pi,phi_obs*180/pi,sqrt(mvar_grd)),colorbar
%contourf(Xgrd,Ygrd,mvar_grd),colorbar
%hold on, %plot([shorelines.Lon],[shorelines.Lat],'k'),
%axis([-50 30 55 82])
%% Sensitivity kernel
%mstd(isnan(mstd))=0;
z0   = 220e3;
meanM = 22e3;
ds   = (maxR-meanM)^2*cos(phi_obs(:))*dlon*dlat;
gszz = zeros(size(Xgrd));
gsxx = zeros(size(Xgrd));
Grr_s = zeros(numel(Xgrd(:)));
Grp_s = zeros(numel(Xgrd(:)));
Gr_s = zeros(numel(Xgrd(:)));
Gphi_s = zeros(numel(Xgrd(:)));
Gxx_s=zeros(ndata,ndata);Gyy_s=Gxx_s;Gzz_s=Gxx_s;Gxy_s=Gxx_s;Gxz_s=Gxx_s;Gyz_s=Gxx_s;
xo = 1e3*Xgrd(:); yo=1e3*Ygrd(:);
rM = maxR-meanM;
[xo3d,yo3d,zo3d]=sph2cart(pi/180*Lon3d(:,:,1),pi/180*Lat3d(:,:,1),Lat3d(:,:,1)*0+rM);
Gzzs = zeros(ndata);
Gxxs=Gzzs;Gyys=Gzzs;Gxys=Gzzs;Gzxs=Gzzs;Gzys=Gzzs;
%
fact = G*drho_s*ds;
%
for i=1:ndata
    %
    rr = sqrt( (xo3d(:)-Xo(i)).^2+(yo3d(:)-Yo(i)).^2+(zo3d(:)-Zo(i)).^2 );
    r2 = rr.*rr;
    r3 = r2.*rr;
    r5 = r3.*r2;
    %
    Gxxs(i,:) = (3*(xo3d(:)-Xo(i)).^2./r5 - 1./r3).*fact;
    Gzzs(i,:) = (3*(zo3d(:)-Zo(i)).^2./r5 - 1./r3).*fact;
    Gyys(i,:) = (3*(yo3d(:)-Yo(i)).^2./r5 - 1./r3).*fact;
    Gxys(i,:) = (3*(xo3d(:)-Xo(i)).*(yo3d(:)-Yo(i))./r5).*fact;
    Gzxs(i,:) = (3*(xo3d(:)-Xo(i)).*(zo3d(:)-Zo(i))./r5).*fact;
    Gzys(i,:) = (3*(yo3d(:)-Yo(i)).*(zo3d(:)-Zo(i))./r5).*fact;
    %
end

%
Cd_zz = (Gzzs*Cms*Gzzs');% variance d2Vdr2
var_zz = diag(Cd_zz)./ cos(phi_obs(:))  ;
%
Cd_xz = (Gzxs*Cms*Gzxs');% variance d2Vdrdphi
var_xz = diag(Cd_xz)./ cos(phi_obs(:))  ;
%
Cd_xx = (Gxxs*Cms*Gxxs');% variance d2Vdrdphi
var_xx = diag(Cd_xx)./ cos(phi_obs(:))  ;
%
Cd_yz = (Gzys*Cms*Gzys');% variance d2Vdrdphi
var_yz = diag(Cd_yz)./ cos(phi_obs(:))  ;
%
Cd_xy = (Gxys*Cms*Gxys');% variance d2Vdrdphi
var_xy = diag(Cd_xy)./ cos(phi_obs(:))  ;
%

Tens_var = [var_xx, var_xy, var_xz, var_yz, var_zz];

figure,
titl = {'Vxx','Vxy','Vxz','Vyz','Vzz'};
for i = 1:size(Tens_var,2)
    cd2d = reshape(Tens_var(:,i),size(phi_obs));
    subplot(3,3,i)
    ax2=axesm('lambertstd','MapLatLimit',[minlat maxlat-2],'MapLonLimit',[minlon maxlon-8],...
        'MLineLocation', 10, 'PLineLocation', 5,...
        'FontSize', 12, 'FontWeight','Bold','LabelFormat','none',...
        'LabelRotation','on','GLineWidth',.5,'GLineStyle','-','GAltitude',4,...
        'MLabelParallel',minlat-.01,'PLabelMeridian',minlon-.1);axis off; framem on; gridm off; mlabel on; plabel on;

    surfm(phi_obs*180/pi,lam_obs*180/pi,cd2d), shading interp
    %c=colorbar(ax2,'Position',[0.7 0.65 0.015 0.15]); ylabel(c, '[E]');
    colorbar
    hold on, plotm([shorelines.Lat],[shorelines.Lon],'k'), shading interp, %colormap(ax2,colbrew(3))
    plotm(plat,plon,'k'),title(titl(i))
    drawnow
end
%
%save DataVariance Cd