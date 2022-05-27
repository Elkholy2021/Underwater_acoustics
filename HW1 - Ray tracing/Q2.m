clear all
clc
close all
%% initialization of inputs
profile_flag=0; %0 if polar region, 1 if the curvy profile
H=3500; %depth
max_time=30; %in seconds
th0=74.17; %angle of sound wave
z0=0; %depth of transducer
x0=0; % x-position of transducer
dt=0.0001; %step size in seconds


%% plotting the profile of sound speed vs depth
tt=max_time;

c0=1450;%initial soundc speed at the ocean surface

z=[0:H];
gradient=1.63e-2;
c_z = @(z) gradient*z+c0;

subplot(1,2,1)
plot(c_z(z),z)

ylim([-1000 H])

set(gca, 'ydir','reverse')
title("Sound speed vs depth")
xlabel('Sound speed')
ylabel('Depth')
grid on
%% Ray trajectory
Z=[];
X=[];
ths=[];
Cs=[];
sign=1;
Xtemp=[];
Ztemp=[];
ctemp=[];
for i =1:tt/dt
    Cs(i)=c0;
    d_v=c0*dt;
    dx=d_v*sind(th0);
    dz=d_v*cosd(th0);
    z=z0+sign*dz;
    x=x0+dx;
    X(i)=x;
    Z(i)=z;
    s=0;
    CCC=c_z(z);

%   th=asind((c_z(z)*sind(th0))/c0);
    th=asind((CCC*sind(th0))/c0);

    if th>90*.99 || z>H*0.9999
        sign=-1; 
    end
    if z < 0.5 
        sign=1;
    end
    ths(i)=th;
    th0=th;
%     c0=c_z(z);
    CCC=c_z(z);

    c0=CCC;
    x0=x;
    z0=z;    
end
subplot(1,2,2)
plot(X/1000,Z)
set(gca, 'ydir','reverse')
title("Ray trajectory")
xlabel('x-direction in kilometer')
ylabel('z-direction in meter -Depth-')
ylim([-1000 H])
grid on