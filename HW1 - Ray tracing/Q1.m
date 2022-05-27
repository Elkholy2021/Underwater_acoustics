clear all
clc
close all
%% initialization of inputs
profile_flag=0; %0 if polar region, 1 if the curvy profile
H=3500; %depth
max_time=30; %in seconds
th0=90; %angle of sound wave
z0=750; %depth of transducer
x0=0; % x-position of transducer
dt=0.0001; %step size in seconds
c0=1450;
zp=[0:1:H];
cp=[];
for i=1:length(zp)
cp(i)=c_z(zp(i),profile_flag);
end
subplot(1,2,1)
plot(cp,zp)
ylim([-1000 H])

set(gca, 'ydir','reverse')
title("Sound speed vs depth")
xlabel('Sound speed')
ylabel('Depth')
grid on

X=[];
Z=[];
sign=1;


for i=1:max_time/dt
    d_v=c0*dt;
    dx=d_v*sind(th0);
    dz=d_v*cosd(th0);
    z=z0+sign*dz;
    x=x0+dx;
    X(i)=x;
    Z(i)=z;
 
    c_new=c_z(z,profile_flag);
%     th=asind((c_new/c0)*sind(th0));
    th=asind((sind(th0)*c_new)/c0);
 

    th0=th;
    c0=c_new;
    x0=x;
    z0=z;
    if th>90*.99 || z>H*0.99
        sign=-1; 
    end
    if z < 0.5 
        sign=1;
    end
    
end
subplot(122)
plot(X/1000,Z)
set(gca, 'ydir','reverse')
title("Ray trajectory")
xlabel('x-direction in kilometer')
ylabel('z-direction in meter -Depth-')
ylim([-1000 H])
grid on
 

function c=c_z(z,flag)
gradient=1.63e-2;
c0=1450;

if flag==0    
    c=gradient*z+c0;
elseif flag==1
    if z <= 700
        c=-gradient*z+c0;
    elseif z>700
        c=gradient*z+c_z(700,flag)-(c0-c_z(700,flag));
    end
end
end
