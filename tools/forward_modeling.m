% forward modeling 
function [Txx_rec,Tyy_rec,Tzz_rec,Tzy_rec,Tzx_rec,Txy_rec,Trr_rec] = forward_modeling(mod_syn,rho_in) 
% ouput in Voight notation
% xx, yy, zz, zy, zx, xy ---> 1, 2, 3, 4, 5, 6 (Cartesian)
% pp, ll, rr, rp, rl, pl ---> 1, 2, 3, 4, 5, 6 (Spherical)
G = 6.67e-11;
h_alt=225e3; 
minR=5971e3;
maxR=6371e3;
dlon=abs(mod_syn.Lon_reg(2)-mod_syn.Lon_reg(1))*pi/180;
dlat=abs(mod_syn.Lat_reg(2)-mod_syn.Lat_reg(1))*pi/180;
dr=abs(mod_syn.r(2)-mod_syn.r(1));
[Lon3d,Lat3d,r3d]=meshgrid(mod_syn.Lon_reg,mod_syn.Lat_reg,mod_syn.r);
[x3d,y3d,z3d]=sph2cart(Lon3d*pi/180,Lat3d*pi/180,r3d);
r_obs = maxR+h_alt;
[Xo,Yo,Zo]=sph2cart(Lon3d(:,:,1)*pi/180,Lat3d(:,:,1)*pi/180,r_obs);
lam = Lon3d(:)*pi/180; lam_obs = Lon3d(:,:,1)*pi/180;
phi = Lat3d(:)*pi/180; phi_obs = Lat3d(:,:,1)*pi/180;
rho_in = rho_in(:); 
nx = length(mod_syn.Lon_reg);
ny = length(mod_syn.Lat_reg);
ndata = nx*ny;
nmod = numel(rho_in);
%% Green functions in Spherical coordinates
%
Gr = zeros(ndata,nmod);Gphi=Gr;Glam=Gr;
Grr=Gr;Grp=Grr;Grl=Grr;Gpp=Grr;Gpl=Grr;Gll=Grr;
L = zeros(ndata,nmod);
%
for i=1:ndata
    cospsi = sin(phi_obs(i))*sin(phi)+cos(phi_obs(i))*cos(phi).*cos(lam_obs(i)-lam);
    %
    L(i,:)  = sqrt(r_obs.^2 + r3d(:).^2 - 2*r3d(:).*r_obs.*cospsi);
    %
    Kphi   = cos(phi_obs(i))*sin(phi)-sin(phi_obs(i))*cos(phi).*cos(lam_obs(i)-lam);
    %
    Gr(i,:) =  r3d(:).^2.*cos(phi)./L(i,:)'.^3 .* (r_obs - r3d(:).*cospsi);
    Gphi(i,:) =  r_obs.*r3d(:).^3 .*cos(phi).* Kphi ./L(i,:)'.^3 ;
    Glam(i,:) =  r_obs*r3d(:).^3.*cos(phi).^2.*cos(phi_obs(i)).*sin(lam_obs(i)-lam) ./L(i,:)'.^3 ;
    %
    Grr(i,:) = r3d(:).^2.*cos(phi)./L(i,:)'.^5 .* (2*r_obs.^2 - r3d(:).^2 - 4*r_obs*r3d(:).*cospsi + 3*r3d(:).^2.*cospsi.^2);
    %
    Grp(i,:) = r3d(:).^3.*Kphi.*cos(phi).* (1 - 3*r_obs.*(r_obs-r3d(:).*cospsi)./L(i,:)'.^2) ./L(i,:)'.^3 ;
    %
    Grl(i,:) = r3d(:).^3.*cos(phi).^2.*cos(phi_obs(i)).*sin(lam_obs(i)-lam) .* ...
        (1 - 3*r_obs.*(r_obs-r3d(:).*cospsi)./L(i,:)'.^2) ./L(i,:)'.^3 ;
    %
    Gpp(i,:)  = r_obs.*r3d(:).^3.*cos(phi) .* ...
        ( 3*r_obs.*r3d(:).*Kphi.^2./L(i,:)'.^2 - cospsi ) ./L(i,:)'.^3 ;
    %
    Gpl(i,:) = r_obs.*r3d(:).^3.*cos(phi).^2.*sin(lam_obs(i)-lam) .*...
        (3*r_obs.*r3d(:).*cos(phi_obs(i)).*Kphi./L(i,:)'.^2 - sin(phi_obs(i))) ./L(i,:)'.^3 ;
    %
    Gll(i,:) = r_obs.*r3d(:).^3.*cos(phi_obs(i)).*cos(phi).^2 .* ...
        ( 3*r_obs.*r3d(:).*cos(phi_obs(i)).*cos(phi).*sin(lam_obs(i)-lam).^2 ./L(i,:)'.^2 - ...
        cos(lam_obs(i)-lam) ) ./L(i,:)'.^3;
end
%% Green functions in Cartesian coordinates
Gzz = zeros(ndata,nmod);
Gxx=Gzz;Gyy=Gzz;Gxy=Gzz;Gzx=Gzz;Gzy=Gzz;
Gx=Gzz;Gy=Gzz;Gz=Gzz;
%
dOm = r3d(:).^2 .* cos(pi/180*Lat3d(:))*dlat*dlon*dr;
dV = dlon*dlat*dr;
%
for i=1:ndata
    %
    rr = sqrt( (x3d(:)-Xo(i)).^2+(y3d(:)-Yo(i)).^2+(z3d(:)-Zo(i)).^2 );
    r2 = rr.*rr;
    r3 = r2.*rr;
        % Gravity vector
    Gx(i,:) =  -(Xo(i)-x3d(:))./r3;
    Gy(i,:) =  -(Yo(i)-y3d(:))./r3;
    Gz(i,:) =  -(Zo(i)-z3d(:))./r3;
    % Gravity gradients tensor Green's functions
    Gxx(i,:) = - 1./r3.*(1 - 3*(Xo(i)-x3d(:)).^2./r2);
    Gyy(i,:) = - 1./r3.*(1 - 3*(Yo(i)-y3d(:)).^2./r2);
    Gzz(i,:) = - 1./r3.*(1 - 3*(Zo(i)-z3d(:)).^2./r2);
    Gzy(i,:) =   1./r3.*(3*(Zo(i)-z3d(:)).*(Yo(i)-y3d(:))./r2);
    Gzx(i,:) =   1./r3.*(3*(Xo(i)-x3d(:)).*(Zo(i)-z3d(:))./r2);
    Gxy(i,:) =   1./r3.*(3*(Xo(i)-x3d(:)).*(Yo(i)-y3d(:))./r2);
end


Txx_rec = G*Gxx*(rho_in.*dOm);
Tzz_rec = G*Gzz*(rho_in.*dOm);
Tyy_rec = G*Gyy*(rho_in.*dOm);
Txy_rec = G*Gxy*(rho_in.*dOm);
Tzx_rec = G*Gzx*(rho_in.*dOm);
Tzy_rec = G*Gzy*(rho_in.*dOm);
%
Vphi_rec= G*Gphi*(rho_in.*dV);
Vlam_rec= G*Glam*(rho_in.*dV);
Vr_rec  = G*Gr*(rho_in.*dV);
Trr_rec = G*Grr*(rho_in.*dV);
Trp_rec = G*Grp*(rho_in.*dV);
Trl_rec = G*Grl*(rho_in.*dV);
Tpp_rec = G*Gpp*(rho_in.*dV);
Tpl_rec = G*Gpl*(rho_in.*dV);
Tll_rec = G*Gll*(rho_in.*dV);
%
Vx_rec = 1./r_obs.*Vphi_rec;
Vy_rec = 1./r_obs./cos(phi_obs(:)).*Vlam_rec;
Vz_rec = Vr_rec;
Vxx_rec = 1./r_obs .* Vr_rec + 1./r_obs.^2 .* Tpp_rec ;
Vyy_rec = 1./r_obs .* Vr_rec - 1./r_obs.^2.*sin(phi_obs(:))./cos(phi_obs(:)).*Vphi_rec + ...
    1./r_obs.^2.*cos(phi_obs(:)).^2 .* Tll_rec;
Vxy_rec = 1./r_obs.^2./cos(phi_obs(:)) .*Tpl_rec + ...
    sin(phi_obs(:))./r_obs.^2./cos(phi_obs(:)).^2 .*Vlam_rec;
Vxz_rec = -1./r_obs.^2 .* Vphi_rec + 1./r_obs.* Trp_rec;
Vyz_rec = -1./r_obs.^2.*cos(phi_obs(:)).*Vlam_rec + 1./r_obs.*cos(phi_obs(:)) .* Trl_rec;
Vzz_rec = Trr_rec;