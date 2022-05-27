clear all
clc
close all
%HW1 - Q3
%% initialization of inputs
profile_flag=0; %0 if polar region, 1 if the curvy profile
H=3500; %depth
max_time=30; %in seconds
th0=90-89.9999999; %angle of sound wave
z0=0; %depth of transducer
x0=0; % x-position of transducer
dt=0.001; %step size in seconds
velocity_at_surface=1450;

%% plotting the profile of sound speed vs depth
tt=max_time;
profile=load('profile.mat'); %loading the curvy profile just in case

if profile_flag==0
    c0=velocity_at_surface;%initial soundc speed at the ocean surface
end
if profile_flag==1
    c0=profile.profile(1,1);%initial soundc speed at the ocean surface
end
z=[0:H];
% gradient=((1520-1450)/4000);
gradient=1.63e-2;
c_z = @(z) gradient*z+c0;

subplot(1,2,1)
if profile_flag==0
    plot(c_z(z),z)
else
    [cd,sf,df]=C_Z(100,H,profile);
    plot(sf,df)
end
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
time_steps=[];
first_hit_surface=0;
for i =1:tt/dt
    time_steps(i)=(i/(tt/dt))*tt;
    Cs(i)=c0;
    d_v=c0*dt;
    dx=d_v*sind(th0);
    dz=d_v*cosd(th0);
    z=z0+sign*dz;
    x=x0+dx;
    X(i)=x;
    Z(i)=z;
    s=0;
    if profile_flag==0
        CCC=c_z(z);
    else
        cd=C_Z2(z,H,profile);
        CCC=cd;
    end
%   th=asind((c_z(z)*sind(th0))/c0);
    th=asind((CCC*sind(th0))/c0);

    if th>90*.99 || z>H*0.99
        sign=-1; 
    end
    if z < 0.5 
        if first_hit_surface==0
            first_hit_surface=1;
            time_first_hit_surface=(i/(tt/dt))*tt
        end
        sign=1;
       
    end
    ths(i)=th;
    th0=th;
%     c0=c_z(z);
    if profile_flag==0
        CCC=c_z(z);
    else
        cd=C_Z2(z,H,profile);
        CCC=cd;
    end
    c0=CCC;
    x0=x;
    z0=z;    
end
disp("Depth if sound velocity is uniform is:")
vpa(velocity_at_surface*(time_first_hit_surface/2))

%subplot(1,2,2)
figure
plot(X*100,Z)
set(gca, 'ydir','reverse')
title("Ray trajectory")
xlabel('x-direction in centimeter')
ylabel('z-direction in meter -Depth-')
ylim([-1000 H])
grid on