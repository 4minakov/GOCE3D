%Conversion of gravity gradients in LNOF to ECEF reference
function [gxx,gyy,gzz,gxy,gzx,gzy]=lnof2ecef(gww,gnn,grr,gwn,grw,grn,LATo,LONo)
gzz = LATo*0; gxx=gzz; gzx = gzz; gyy = gzz;
gzy = gzz; gxy = gzz;
[nxo,nyo]=size(LATo);
for i = 1:length(LATo(:))
    az1=LONo(i)*pi/180; el1=LATo(i)*pi/180;
    a=[-sin(az1), cos(az1), 0;...
        -sin(el1)*cos(az1),-sin(el1)*sin(az1),cos(el1);...
        cos(el1)*cos(az1),cos(el1)*sin(az1),sin(el1)];
    v1 = [gww(i),-gwn(i),-grw(i);...
        -gwn(i),gnn(i),grn(i);...
        -grw(i),grn(i),grr(i)];
    vs1=a'*v1*a;
    gxx(i) =  vs1(1,1); gyy(i) =  vs1(2,2);
    gxy(i) =  vs1(1,2); gzy(i) =  vs1(2,3);
    gzx(i) =  vs1(1,3); gzz(i) =  vs1(3,3);
  
end
gzz = reshape(gzz,nxo,nyo); gxx = reshape(gxx,nxo,nyo);
gzx = reshape(gzx,nxo,nyo); gzy = reshape(gzy,nxo,nyo);
gyy = reshape(gyy,nxo,nyo); gxy = reshape(gxy,nxo,nyo);