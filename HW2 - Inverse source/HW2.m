clear all
clc
close all
load("received_signals.mat");
t=linspace(-0.6252,2.5-0.6252,3001);
t2=linspace(-0.6252,2.5*2-0.6252,3000*2);
c=1500;
sum=zeros(1,3001);
reversed_signals=[];
h=150;
resolution = 10; 
yRange = (0 : resolution : h); 
xRange = (0 : resolution : 1500); 
hr = 15 : 15 : h - 15; 
ns = 8;

%% Flipping signals
for s = 1:9
subplot(330+s)
plot(t,green(s,:))
hold on
reversed_signals(s,:)=flip(green(s,:));
hold on
% plot (t,reversed_signal)
end
%% Backpropagation
grid=[];
for xn =1: length(xRange)
    for zn =1: length(yRange)
        xs=xRange(xn);
        zs=yRange(zn);
        xr=1200;
        received_signal=zeros(1,6000);
        for i = 1:9
            zr=hr(i);
            signal=reversed_signals(i,:);
            received_signal = apply_green_function(signal,c,h,ns,xr,zr,xs,zs)+received_signal;
        end
        received_signal=normalize(received_signal);
        grid(zn,xn)=max(received_signal);
    end
end
figure
image(grid,'CDataMapping','scaled') 
xticklabels((xticks*resolution)-resolution+xRange(1)) 
yticklabels((yticks*resolution)-resolution+yRange(1))
colorbar
[M,I]=max(grid);
[M2,I2]=max(M);
depth=I(I2);
x_distance=I2;
disp("Source depth is")
yRange(depth)
disp("Source x-distance is")
xr-xRange(x_distance)
title("Source is "+compose("%d",round(yRange(depth)))+"m from surface and "+compose("%d",round(xr-xRange(x_distance)))+"m from the recievers")
