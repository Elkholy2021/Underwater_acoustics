clear variables;
clc
close all;

%% Inputs
h = -3500; %height to sea bottom
grad_inv_h = -700; %height at which speed gradient will be inverted
zs = -700;
th = 89;
gradient = 1.63e-2; %speed gradient
dt = 0.0001;
max_time = 30;
c0 = 1450; %base speed for c
h_Range = -(0:1:abs(h));
speedProfile = arrayfun(@(h) getC(h,gradient,c0, grad_inv_h),h_Range);
time = (0:1:(max_time/dt));

%% Script

% points = getPoints(zs, th ,dt, gradient, max_time, c0, h, grad_inv_h);
% points2 = getPoints(zs, th + 7 ,dt, gradient, max_time, c0, h, grad_inv_h);
% points3 = getPoints(zs, th + 15 ,dt, gradient, max_time, c0, h, grad_inv_h);
% subplot(1,4,2:4);
% plot(points(:,1), points(:,2))
% title("Ray Trayectory x,z")
% subplot(1,4,1);
% plot(arrayfun(@(h) getC(h,gradient,c0, grad_inv_h),h_Range), h_Range) 
% title("Speed Profile")


%% Functions
% function tracingPoints = getPoints(zs, thOrigin, dt, grad, totalT, cStart, h, grad_inv_h)
x = 0;
z = zs;
t = 0;
tracingPoints = zeros(max_time/dt, 4);
th =  th ;
C = getC(z,gradient,c0, grad_inv_h);
index = 1;

while t < max_time
    % compute step
    [dz, dx]= stepDt(th, dt, C);

    % Log values
    tracingPoints(index,1) = x;
    tracingPoints(index,2) = z;
    tracingPoints(index,3) = C;
    tracingPoints(index,4) = th;

    % Check reflection
    if z + dz > 0 || z + dz < h
        th = 180 - th;
        x = x + dx;
        t = t + dt;
        index = index + 1;
        continue
    end

    % Get new C and th
    newC = getC((z + dz), gradient, c0, grad_inv_h);
    newTheta =  getTheta(newC, C, th);

    % Update values for next iteration
    C = newC;
    th = newTheta;
    t = t + dt;
    x = x + dx;
    z = z + dz;
    index = index + 1;
end
% end
subplot(1,4,2:4);
plot(tracingPoints(:,1), tracingPoints(:,2))
title("Ray Trayectory x,z")
subplot(1,4,1);
plot(arrayfun(@(h) getC(h,gradient,c0, grad_inv_h),h_Range), h_Range) 
title("Speed Profile")

function newTheta = getTheta(newC, C, th)

if abs(th - 90) < 0.01 % check if ray is aproaching 90 from below or above
    th = th - sign(th - 90)*0.02; % bump ray depending on the direction it was coming from
end

% asin((C_i / C_i-1) * sin(th)) function is symetrical and mirrored at
% 90 degrees, angles beyond 90 degrees have to be translated to the
% opposite plane in order for the output to match the input
if th > 90
    newTheta = 180 - asind((newC / C) * sind(th));
else
    newTheta = asind((newC / C) * sind(th));
     
end
end

function [dz, dx] = stepDt(th, dt, c)
cdt = c * dt;
dz = -cdt*cosd(th);
dx = cdt*sind(th);
end

function c = getC(z, grad, base, h_inv)
if z >= h_inv %this keeps the gradient function  continuous
    c = base + (-grad * abs(z)) + (grad * abs(h_inv)) - (-grad * abs(h_inv));
else
    c = base + (grad * abs(z));
end
end