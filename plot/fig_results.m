% fig_results
%
% Scripts to calculate and plot figures as in manuscript by
%
% Minakov, A., & Gaina, C. (2021).
% Probabilistic linear inversion of satellite gravity gradient data applied
% to the northeast Atlantic. Journal of Geophysical Research: Solid Earth,
% 126, e2021JB021854. https://doi.org/10.1029/2021JB021854
% 
% Last modified by alexamin@uio.no, 26/11/2021
%
% version v1.1
% 
% Dataset in ..data/GOCE_NEATLANTIC is structure containing the full model
% 
%      Cm: [6670×6670 double] posterior model covariance matrix
%       m: [29×23×10 double] mean denstity perturbation model
%      Cd: [667×667 double] data covariance matrix
%       d: [29×23 double] data vector (Trr)
%       r: [1×10 double] distance
%     lat: [29×1 double] latitute
%     lon: [23×1 double] longitude

% for plotting some files GMT software needs to to installed
% and paths added 
%
clear all 
close all
%
% Loading Datasets

LoadingDatasets

% Processing

ProcessingResults


% Plotting

fig1

fig2

fig3

fig4

fig5
%
fig6

fig7

fig8

fig9gmt

fig10gmt

fig11gmt

fig12

fig13gmt

fig14gmt
%%
fig15

fig16

fig17

fig18gmt

fig19

fig20

fig21

fig22

fig23

fig24

