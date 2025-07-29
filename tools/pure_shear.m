function Tp = pure_shear(betaf,Age,Z)
%Plate cooling model
Tm = 1;
L = 90e3;
k = 1e-6;
tau = L^2/pi^2/k;
Z = Z*1e3;
%ind = find(Z<=L);
%Tp = zeros(size(Z))+Tm;
if Z <= L
    %Temperature
    const = (2*Tm)/pi;
    N=10;
    sumtotal = 0;
    for n=1:N
        r = (betaf./(n*pi)).*sin((n*pi)./betaf);
        sum = (((-1)^(n+1))/n)*r.*(exp(((-(n^2).*Age)./tau))).*sin(n*pi.*(1-Z/L));
        sumtotal = sum + sumtotal;
    end
    Tp = Tm*Z/L + const*sumtotal;
else
    Tp = Tm+zeros(size(Age));
end
%%
