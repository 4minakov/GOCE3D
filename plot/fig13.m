%figure 13 CSEM model

if ~isfile('../data/csem-north-atlantic-2019.12.01.nc')
    disp('Downloading tomography model ..')
    websave('../data/csem-north-atlantic-2019.12.01.nc','http://ds.iris.edu/files/products/emc/emc-files/csem-north-atlantic-2019.12.01.nc')
end
%
%addpath('C:\Program Files\gmt6\bin')

tomo_info = ncinfo('csem-north-atlantic-2019.12.01.nc');
tomoNA.depth = ncread(tomo_info.Filename,'depth');
tomoNA.latitude = ncread(tomo_info.Filename,'latitude');
tomoNA.longitude = ncread(tomo_info.Filename,'longitude');
tomoNA.vsv = ncread(tomo_info.Filename,'vsv');
tomoNA.vsh = ncread(tomo_info.Filename,'vsh');
v3d = double(tomoNA.vsv); v3d=v3d-mean(mean(v3d,1),2);
sz = size(v3d);
T1 = v3d;
R1 = double(tomoNA.depth);

ET=sum(sum(sum(T1(:,:,7:40))))/prod(sz(:));     %mean
T1=T1-ET; %demean
VT=sum(sum(sum(T1(:,:,7:40).^2)))/prod(sz(:));  %dispersion

%% radial covariance estimation
nbin = 30;
Vr=zeros(1,nbin);           % variogram in radial direction
dd = linspace(0,400,nbin);
for id = 2:nbin-1
    disp(id)
    m=0;
    sm=0;
    for  i=7:40
        for j=7:40
            dst = abs(R1(i)-R1(j));
            if (dst>dd(id-1)) && (dst<dd(id+1)) && (i~=j)
                tt1 = T1(:,:,i);
                tt2 = T1(:,:,j);
                sm = sm + sum((tt1(:)-tt2(:)).^2);
                m  =  m + sz(1)*sz(2);
            end
        end
    end
    Vr(id)=0.5*sm/m;     
end
Cr = VT - Vr;         
%fileout = ['../data/','CorrRad_',num2str(nbin),'bin'];
%save(fileout,'Cr','dd')
%% angular covariance estimation
nbin = 30;
Vt=zeros(1,nbin);           % variogram in angular direction
dda = linspace(0,30*pi/180,nbin);
%llat = vrick.Lat3d(:,:,1);
%llon = vrick.Lon3d(:,:,1);
[llat,llon]=meshgrid(tomoNA.latitude,tomoNA.longitude);
iz  = 15;
Ct = 0;
step = 8;
for id = 2:nbin-1
    disp(id)
    m=0;
    sm=0;
    for i1=1:step:sz(1)
        for j1=3:step:sz(2)
            for i2=2:step:sz(1)
                for j2=1:step:sz(2)
                    if (i1~=i2) || (j1~=j2)
                        dst = pi/180*distance(llat(i1,j1),llon(i1,j1),llat(i2,j2),llon(i2,j2));
                        tt1 = T1(i1,j1,iz);tt1=tt1(:);
                        tt2 = T1(i2,j2,iz);tt2=tt2(:);
                        if dst>dda(id-1) && dst<=dda(id+1)
                            sm = sm + sum((tt1-tt2).^2);
                            m  =  m + 1;
                        end
                    end
                end
            end
        end
    end
    if m>0
        Vt(id)=0.5*sm/m;     
    end
end
Vt = Vt/max(Vt);
Ct = Ct + (1 - Vt);           % angular covariance

%fileout = ['../data/','CorrAng_',num2str(ceil(R1(iz))),'km'];
%save(fileout,'Ct','dda')
%% radial covarinace theoretical
xh = linspace(0.1,400,100);
df = 1./(xh(2)-xh(1));
f = df*[0:length(xh)/2];
n = 2.5;
alpha_b = 1/0.1e3;
alpha_e = 1/.2e3;
w0 = 2*pi/.6e3;
tau   = xh;
C_e = 1;
C_b = 1;
B_exp = C_e*exp(-alpha_e*tau).*cos(w0*tau);
%B_bes = C_b*(alpha_b*tau).^(-(n-2)/2).*besselj((n-2)/2,alpha_b*tau);
%B_bes = B_bes/B_bes(1);

% angular covariance theoretical
th = linspace(0,0.5,100);
C_a1 = 0; 
C_a2 = 0;
ro1 = .9;
ro2 = .99;
for n = 0:51
    Pl = legendre(n,cos(th)); Pl = Pl(1,:);
    C_a1 = C_a1 + (ro1^n)*Pl;
    C_a2 = C_a2 + (ro2^n)*Pl;
end

%%
figure
subplot(121),
p1=plot(dd,Cr/Cr(1),'-o','LineWidth',1); hold on, 
plot(xh,B_exp,'-k','LineWidth',1)
xlabel('Distance (km)'), ylabel('Correlation'), xlim([0 380]), 
tt = title(['(a)      Radial ']); tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left'; set(gca,'FontSize',12,'LineWidth',1)
subplot(122),
p1=plot(dda,Ct,'o-','LineWidth',1); hold on, 
plot(th,C_a1/C_a1(1),'k','LineWidth',1),plot(th,C_a2/C_a2(1),'k','LineWidth',1),
xlabel('Radian'), xlim([0 0.5])
tt = title(['(b)     Angular ']); tt.Units='Normalized'; tt.Position(1)=0;
tt.HorizontalAlignment='left'; set(gca,'FontSize',12,'LineWidth',1)
text(0.3,0.4,'\rho = 0.9','FontSize',11), text(0.05,0.1,'\rho = 0.99','FontSize',11), 
set(gcf,'units','normalized','outerposition',[0.1 0.1 0.65 .6])
print('-depsc','-r400','../fig/fig13')

%gmt('psconvert','../fig/fig13.eps -Tf -P -A ')
%open('../fig/fig13.pdf')

