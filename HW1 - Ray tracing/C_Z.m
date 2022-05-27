function [cd,sf,df]=C_Z(zd,H,profile)
%     zd=0.4857;
%     H=15000;
%     profile=load('profile.mat');
    
    profile=profile.profile;
    n=500;
    speed=profile(:,1);
    depth=profile(:,2);
    s=linspace(speed(1),speed(2),n);
    d=linspace(depth(1),depth(2),n);

    for i=2:length(speed)-1
        sn=linspace(speed(i),speed(i+1),n);
        s=horzcat(s,sn);
        dn=linspace(depth(i),depth(i+1),n);
        d=horzcat(d,dn);

    end
    m=85;
    z=linspace(d(length(d)),H,20000);
    C=@(z) (1/m)*(z-d(length(d)))+s(length(s));
    CC=C(z);
    s=horzcat(s,CC);
    d=horzcat(d,z);
    sf=s;
    df=d;

    for i=1:length(depth)
        if depth(i) > zd
            id=i;
            zi=depth(i);
            break
        end
    end
    if depth(i)<  zd
        id=i;
        zi=depth(i);
    end
if id >1
        for i=1:n
            if d(n*(id-2)+i)>zd
                id2=i;
                zi2=d(n*(id-2)+i);
                cd=(s(n*(id-2)+i)+s(n*(id-2)+i-1))/2;
                break
            end
        end
        if depth(length(depth))< zd
            cd=C(zd);
        end
end
if id==1
    for i=1:n
    if d(n*(id)+i)>zd
        id2=i;
        zi2=d(n*(id)+i);
        cd=(s(n*(id)+i)+s(n*(id)+i-1))/2;
        break
    end
    end
    if depth(length(depth))< zd
    cd=C(zd);
    end
end
end
