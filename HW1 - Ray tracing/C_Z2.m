function [cd]=C_Z(zd,H,profile)
% clear all
% clc
% close all
% zd=500;
% H=15000;
% profile=load('profile.mat');


    profile=profile.profile;
    n=500;
    m=85;
    speed=profile(:,1);
    depth=profile(:,2);
    if zd> depth(length(depth))
        C=@(z) (1/m)*(z-depth(length(depth)))+speed(length(speed));
        CC=C(zd);
    end

    if zd <depth(length(depth))
        for i=1:length(depth)
            if depth(i) > zd
                id=i;
                zi=depth(i);
                break
            end
        end
        if id >1
            m=(depth(id)-depth(id-1))/(speed(id)-speed(id-1));
            C=@(z) (1/m)*(z-depth(id-1))+speed(id-1);
            CC=C(zd);

        end

        if id==1
            CC=speed(id);
        end

    end
    cd=CC;
end