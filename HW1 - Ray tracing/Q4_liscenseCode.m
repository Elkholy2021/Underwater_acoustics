clear all
clc
close all
%% Q1
% xo=0;
% zo=750;
% theta0=-0.001;
% tt=30;
% zz=linspace(0,3500,100000);
% cc=1450+zz.*1.63*10e-2;
% plotflag=1;
%% Q2
% xo=0;
% zo=0;
% theta0=90-74.25;
% tt=30;
% zz=linspace(0,3500,100000);
% cc=1450+zz.*1.63e-2;
% plotflag=1;
%% Q4
xo=0;
zo=700;
theta0=2;
tt=30;
zz1=linspace(0,700,100000/5);
zz2=linspace(701,3500,4*100000/5);
zz=cat(2,zz1,zz2);
cc1=1450+zz1.*-1.63*10e-2;
cc2=cc1(length(cc1))-(1450-cc1(length(cc1)))+zz2.*1.63*10e-2;
cc=cat(2,cc1,cc2);
plotflag=1;


%%
[xxf, zzf, ttf, ddf] = raytrace(xo,zo, theta0,tt,zz,cc,plotflag)
