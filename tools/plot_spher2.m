%% Plotting
function plot_spher2(T,Lat_reg,Lon_reg,r)
%figure,
%testing 
% [TH,PH,R0] = ndgrid(pi/2-theta,phi,R);
% [XX,YY,ZZ] = sph2cart(PH,TH,R0);
% T = exp(-( (XX-0.7).^2+(YY-0.7).^2+(ZZ).^2 )/0.5^2);
%%
theta = pi/2-Lat_reg*pi/180;
phi = Lon_reg*pi/180;
R = 2*r-6371e3;
R1 = 6371e3-2*300e3;
R2 = 6371e3;
%T = m_out3d;
NR = length(r);
Ntheta=length(Lat_reg);
Nphi=length(Lon_reg);
%plot values for R_top
[TH,PH] = ndgrid(pi/2-theta,phi);
T2d = squeeze(T(:,:,NR-1));
R0 = R2*ones(size(TH));
%ind = find(PH>0&PH<2*pi/3&TH>0);
T2d1 = T2d; %T2d1(ind)=NaN;
[X,Y,Z] = sph2cart(PH,TH,R0);
surf(X,Y,Z,T2d1), axis equal tight
shading flat, view(148,28)
hold on,
%%
%plot values for R_bottom
T2d = squeeze(T(:,:,2));
R0 = R1*ones(size(TH));
[X,Y,Z] = sph2cart(PH,TH,R0);
surf(X,Y,Z,T2d),shading flat,
%%
% plot values for min theta
[PH,R0] = ndgrid(phi,R);
T2d = squeeze(T(theta==theta(1),:,:));
[X,Y,Z] = sph2cart(PH,PH*0+pi/2-theta(1),R0);
surf(X,Y,Z,T2d),shading flat
%%
% plot values for min phi
[TH,R0] = ndgrid(pi/2-theta,R);
T2d = squeeze(T(:,1,:));
[X,Y,Z] = sph2cart(TH*0+phi(1),TH,R0);
hold on
surf(X,Y,Z,T2d),shading flat
%%
% plot values for max theta
[PH,R0] = ndgrid(phi,R);
T2d = squeeze(T(theta==theta(Ntheta),:,:));
[X,Y,Z] = sph2cart(PH,PH*0+pi/2-theta(Ntheta),R0);
surf(X,Y,Z,T2d),shading flat
%%
% plot values for max phi
[TH,R0] = ndgrid(pi/2-theta,R);
T2d = squeeze(T(:,Nphi,:));
[X,Y,Z] = sph2cart(TH*0+phi(Nphi),TH,R0);
hold on
surf(X,Y,Z,T2d),shading flat
%%
camlight, material dull, lighting Gouraud, shading interp
axis off
%%
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 4]),
%print('T_vertical','-dpng','-r600')