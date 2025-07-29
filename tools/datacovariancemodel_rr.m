%% data covariance matrix
G = 6.67e-11;
phi0 =  90;%true latitude
lam0 = -0;%zero longitude
mkrig_var = load('moho_var.out');
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
xo = 1e3*Xgrd(:); yo=1e3*Ygrd(:);
rM = maxR-meanM;
for i = 1:length(phi_obs(:))   
    cospsi = sin(phi_obs(i))*sin(phi_obs(:))+...
        cos(phi_obs(i))*cos(phi_obs(:)).*...
        cos(lam_obs(i)-lam_obs(:));
    %
    L = sqrt(r_obs.^2 + rM.^2 - 2*rM.*r_obs.*cospsi);
    %
    %Gr   =  -G * drho_s.* (r_obs - rM.*cospsi) ./L.^3 .* ds ;%gravity (dVdr)
    %
    Grr_s(i,:) = G * drho_s ./L.^5 .* (2*r_obs.^2 - rM.^2 - 4*r_obs*rM.*cospsi + 3*rM.^2.*cospsi.^2) .* ds;
    %Grr = G * drho_s * (- 1./L.^3 + 3*(r_obs - rM*cospsi).^2./L.^5) .* ds ;
    
end

%Cd_r = (Gr_s*Cms*Gr_s');
%var_r = diag(Cd_r)./ cos(phi_obs(:))  ;
%
Cd_rr = (Grr_s*Cms*Grr_s');% variance d2Vdr2
var_rr = diag(Cd_rr)./ cos(phi_obs(:))  ;
%save DataVariance Cd_rr