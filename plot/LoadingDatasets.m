if ~isfolder('../tools/SHbundle-master')
    disp('Downloading SHbundle m-files ..')
    websave('../tools/SHbundle-master.zip','https://www.gis.uni-stuttgart.de/dokumente/SHbundle-master.zip')
    unzip('../tools/SHbundle-master.zip','../tools/SHbundle-master')
    system('del ..\tools\*.zip')
end
if ~isfolder('../tools/uberall-master')
    disp('Downloading SHbundle m-files ..')
    websave('../tools/uberall-master.zip','https://www.gis.uni-stuttgart.de/dokumente/uberall-master.zip')
    unzip('../tools/uberall-master.zip','../tools/uberall-master')
    system('del ..\tools\*.zip')
end
if ~isfolder('../tools/ScientificColourMaps7')
    disp('Downloading colormaps ..')
    websave('../tools/ScientificColourMaps7.zip','https://zenodo.org/record/5501399/files/ScientificColourMaps7.zip')
    unzip('../tools/ScientificColourMaps7.zip','../tools/ScientificColourMaps7')
    system('del ..\tools\*.zip')
end
if ~isfile('../data/XGM2016.gfc')
    disp('Downloading XGM gravity model ..')
    websave('../data/XGM2016.zip','https://datapub.gfz-potsdam.de/download/10.5880.ICGEM.2017.003/XGM2016.zip')
    unzip('../data/XGM2016.zip','../data')
    system('del ..\data\*.zip')
end
%
addpath(genpath('../data'))
addpath(genpath('../tools'))

load XGM
load topography
load ice
load sediments
load Moho_combined

load OceanAge
age1=age;age1(age<0)=max(age(:));

load model_PS_new90

Vp_data = load('szwillus_vp.txt');
% MLon=reshape(Vp_data(:,1),360,180);
% MLat=reshape(Vp_data(:,1),360,180);
M_data = load('szwillus_moho.txt');
M = reshape(M_data(:,3),360,180);
MLat = reshape(M_data(:,2),360,180);
MLon = reshape(M_data(:,1),360,180);
Vp = reshape(Vp_data(:,3),360,180);

load GSHHS_i
load roma

load figData 
load rwbcmap
load bwcmap