%birds a la vicsek et al PRL 1995
% not designed to be efficient or fast!
clear all
close all
lbox=25.
%vs=0.03;
vs=0.5;
dt=1.;
%interaction radius
r=12.5;
%r=0.2;
mindex=0;
%random fluctuation in angle per step (rads)
%order
%eta=(2.*pi).*.005
eta=1;
% near critical
%
%eta=(2.*pi).*.3
%disorder
%eta=20
px=NaN(8,5000);
py=px;
vx=px;
vy=px;
mangui=px;
nbird=8;
birdl=[1:nbird];
onebd=ones(1,nbird)';
%bird arrays position and angle
figure
axis([0 lbox 0 lbox])
axis('square')
hold on

%IC random

xb=rand(nbird,1).*lbox;
yb=rand(nbird,1).*lbox;

ang=pi.*2.*rand(nbird,1);
nsteps=1;
vxb=vs.*cos(ang);
vyb=vs.*sin(ang);
vx(:,nsteps)=vxb;
vy(:,nsteps)=vyb;
%outer loop

px(:,nsteps)=xb;
py(:,nsteps)=yb;
for nsteps=2:5000;
    nsteps;
   
    
xb=xb+vxb.*dt;
yb=yb+vyb.*dt;
px(:,nsteps)=xb;
py(:,nsteps)=yb;
for bird1=1:nbird;
%periodic
if(xb(bird1)<0);xb(bird1)=xb(bird1)+lbox; end
if (yb(bird1)<0);yb(bird1)=yb(bird1)+lbox;end
if (xb(bird1)>lbox);xb(bird1)=xb(bird1)-lbox;end
if (yb(bird1)>lbox);yb(bird1)=yb(bird1)-lbox;end

%find mean angle of neigbours (include bird1)

    sep(1:nbird)=sqrt((onebd.*xb(bird1)-xb(1:nbird)).^2+...
        (onebd.*yb(bird1)-yb(1:nbird)).^2);
    xseq=[1 2 3 4 5 6 7 8];
    nearang=xseq(sep<r);
    vxtemp=mean(vx(nearang,nsteps-1));
    vytemp=mean(vy(nearang,nsteps-1));
    vtemp=[vxtemp,vytemp];
    if vytemp>=0
    costheta(bird1)=acos(vxtemp/norm(vtemp));
    else if vytemp<0
    costheta(bird1)=2*pi-acos(vxtemp/norm(vtemp));       
%     vx(bird1,nsteps)=vxb;
%     vy(bird1,nsteps)=vyb;
    %mang(bird1)=mean(nearang);
        end
    end
    
end
    ang=costheta'+eta.*(rand(nbird,1)-0.5);
    %ang=costheta';
    if ang<0 ang=ang+2*pi;end
    if ang>2*pi ang=ang-2*pi;end

vxb=vs.*cos(ang);
vyb=vs.*sin(ang);
vx(:,nsteps)=vxb;
vy(:,nsteps)=vyb;
%mangui(:,nsteps)=mang;
% teh order parameter

%vtot(nsteps)=abs(mean(vxb).^2+mean(vyb).^2);
%dist(nsteps)=mean(sep);

%plot
cla

set(gcf,'doublebuffer','on')

plot(xb,yb,'b.')
%quiver(xb,yb,vxb,vyb,'b')

drawnow

    
end
vx1=vx.^2;
vy1=vy.^2;
d=vx1+vy1;
d=sqrt(d);
%vx2 is velocity unitization
vx2=vx./d;
vy2=vy./d;
%Calculate the degree of polarization
vx3=sum(vx2)/8;
vy3=sum(vy2)/8;
%p is degree of polarization
p=sqrt(vx3.^2+vy3.^2);
vgx=vx3./p;
vgy=vy3./p;
h=figure;
plot(p)
% %% jisuan jiaosulv
% for i=2:5000
%     vcos(i)=vgx(i)*vgx(i-1)+vgy(i)*vgy(i-1);
% end
% ac=acos(vcos)/pi*180;
% h=figure;
% hist(ac,500);
% save('/Users/zq/Desktop/file/20150902/viseck_data_0.5.mat');