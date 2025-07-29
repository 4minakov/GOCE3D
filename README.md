# GOCE3D

MATLAB scripts to calculate and plot figures as in manuscript by

 Minakov, A., & Gaina, C. (2021).
 Probabilistic linear inversion of satellite gravity gradient data applied
 to the northeast Atlantic. Journal of Geophysical Research: Solid Earth,
 126, e2021JB021854. https://doi.org/10.1029/2021JB021854
 
 Last modified by alexamin@uio.no, 26/11/2021

 version v1.1
 

 Contents of arhcive
 /data  contains requiried and generated datasets 
 /fig   folder for output figures 
 /plot  scripts to produce figures 
 /tools additional matlab tools and routines

 Dataset in ..data/GOCE_NEATLANTIC is structure containing the full model
 
      Cm: [6670×6670 double] posterior model covariance matrix
       m: [29×23×10 double] mean denstity perturbation model
      Cd: [667×667 double] data covariance matrix
       d: [29×23 double] data vector (Trr)
       r: [1×10 double] distance
     lat: [29×1 double] latitute
     lon: [23×1 double] longitude

 Run  /plot/fig_results.m to produce all figures 

 Some scripts require GMT (Wessel et al. 2019) and SHBUNDLE (Sneeuw et al. 2018) software to be installed

 and corresponding folders must be added to the matlab search path.

 Also ScientificColorMaps7 by F. Crameri (2021) maybe required and have been included in the archive.
