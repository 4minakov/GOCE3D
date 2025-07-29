function [x,y] = pstereo(phi,lam,phi0,lam0)
R = 6371;
phi = phi*pi/180; lam = lam*pi/180;
phi0 = phi0*pi/180;
lam0 = lam0*pi/180;
k = 2*R./(1+ sin(phi0).*sin(phi)+cos(phi0).*cos(phi).*cos(lam-lam0));
x = k.*cos(phi).*sin(lam-lam0);
y = k.*(cos(phi0).*sin(phi)-sin(phi0).*cos(phi).*cos(lam-lam0));