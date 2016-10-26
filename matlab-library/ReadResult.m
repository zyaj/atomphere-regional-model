%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ReadResult.m
% usage: read the results from the output 
% of regional model for cee263B 
% save as result.mat which can be used for
% further discussion
% Yun Zhang 05/06/2015
% @Stanford
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% result directory
datadir='../results';

NLAT=40;
NLONG=40;
NVERT=15;
lat_0=-1;
long_0=-1;
dlat=0.05;
dlong=0.05;
Pa_top=250;
for i=1:NLAT
    for j=1:NLONG
        lat(i,j)=lat_0+(i-1)*dlat;
        long(i,j)=long_0+(j-1)*dlong;
    end
end

latuface=zeros(NLAT,NLONG+1);
longuface=zeros(NLAT,NLONG+1);
latuface(:,1:NLONG)=lat;
latuface(:,NLONG+1)=lat(:,1);
longuface(:,1:NLONG)=long-0.5*dlong;
longuface(:,NLONG+1)=long(:,NLONG)+0.5*dlong;

latvface=zeros(NLAT+1,NLONG);
longvface=zeros(NLAT+1,NLONG);
longvface(1:NLAT,:)=long;
longvface(NLAT+1,:)=long(NLAT,:);
latvface(1:NLAT,:)=lat-0.5*dlong;
latvface(NLAT+1,:)=lat(NLAT,:)+0.5*dlong;

% load in all results
filename=[datadir,'/pi.txt'];
pi_tmp=load(filename);

filename=[datadir,'/Pa.txt'];
Pa_tmp=load(filename);

filename=[datadir,'/geopot.txt'];
geopot_tmp=load(filename);

filename=[datadir,'/K_t.txt'];
K_t_tmp=load(filename);

filename=[datadir,'/nu.txt'];
nu_t_tmp=load(filename);

filename=[datadir,'/PVT.txt'];
PVT_tmp=load(filename);

filename=[datadir,'/u.txt'];
u_tmp=load(filename);

filename=[datadir,'/v.txt'];
v_tmp=load(filename);

filename=[datadir,'/w.txt'];
w_tmp=load(filename);

filename=[datadir,'/qv.txt'];
qv_tmp=load(filename);

filename=[datadir,'/temp.txt'];
temp_tmp=load(filename);

filename=[datadir,'/rhoa.txt'];
rhoa_tmp=load(filename);

filename=[datadir,'/gas.txt'];
gas_tmp=load(filename);

% get Nt from output
Nt=length(pi_tmp')/NLAT;

% seperate different time step
for i=1:Nt
  base1=(i-1)*NLAT;
  base2=(i-1)*NLAT*NLONG;
  base3=(i-1)*(NLAT+1)*NLONG;
  base4=(i-1)*(NLONG+1)*NLAT;
  pi{i}=pi_tmp((base1+1):(base1+NLAT),:);
  Pa{i}=Pa_tmp((base2+1):(base2+NLAT*NLONG),:);
  geopot{i}=geopot_tmp((base2+1):(base2+NLAT*NLONG),:);
  K_t{i}=K_t_tmp((base2+1):(base2+NLAT*NLONG),:);
  nu_t{i}=nu_t_tmp((base2+1):(base2+NLAT*NLONG),:);
  PVT{i}=PVT_tmp((base2+1):(base2+NLAT*NLONG),:);
  u{i}=u_tmp((base4+1):(base4+NLAT*(NLONG+1)),:);
  v{i}=v_tmp((base3+1):(base3+(NLAT+1)*NLONG),:);
  w{i}=w_tmp((base2+1):(base2+NLAT*NLONG),:);
  qv{i}=qv_tmp((base2+1):(base2+NLAT*NLONG),:);
  temp{i}=temp_tmp((base2+1):(base2+NLAT*NLONG),:);
  rhoa{i}=rhoa_tmp((base2+1):(base2+NLAT*NLONG),:);
  gas{i}=gas_tmp((base2+1):(base2+NLAT*NLONG),:);
end

% save results
save('results.mat');

% calculate time accuracy
load('results dt0.625 backward.mat');
aa(1,:)=u{i}(800,:);
load('results dt1.25 backward.mat');
aa(2,:)=u{i}(800,:);
load('results dt2.5 backward.mat');
aa(3,:)=u{i}(800,:);
load('results dt5 backward.mat');
aa(4,:)=u{i}(800,:);
load('results dt10 backward.mat');
aa(5,:)=u{i}(800,:);

load('results dt0.625 matsuno.mat');
bb(1,:)=u{i}(800,:);
load('results dt1.25 matsuno.mat');
bb(2,:)=u{i}(800,:);
load('results dt2.5 matsuno.mat');
bb(3,:)=u{i}(800,:);
load('results dt5 matsuno.mat');
bb(4,:)=u{i}(800,:);
load('results dt10 matsuno.mat');
bb(5,:)=u{i}(800,:);

for i=1:4
  err1(:,i)=aa(:,i)-aa(:,i+1);
  err2(:,i)=bb(:,i)-bb(:,i+1);
end

err11=max(abs(err1));
err22=max(abs(err2));

Nt=11

% plot pi
for i=1:Nt
% figure(1)
% plot2Dresults(pi{i},long,lat,'longitude (degree)','latitude (degree)','column pressure field (hpa)')
% caxis([759.85 759.99])
%figure(2)
 plot3Dsliceresults(log(abs(gas{i})),1,long,lat,'longitude (degree)','latitude (degree)','gas concentration (log(kg/m3))')
 caxis([-8 0])
% figure(3)
% plot3Dsliceresults(v{i},1,longvface,latvface,'longitude (degree)','latitude (degree)','v field (m/s)')
% 
% figure(4)
% plotvelocityquiver(u{i},v{i},7,long,lat,'longitude (degree)','latitude (degree)','velocity vector field')
% axis([min(min(long)),max(max(long)),min(min(lat)),max(max(lat))])
a=num2str(i);
print(gcf,'-r250','-dpng',a)
clf
end

% plot Pa_surf
plot2Dresults(pi{2}+Pa_top,lat,long,'latitude (degree)','longitude (degree)','surface pressure field (hpa)')

% Plot Pa_c
plot3Dsliceresults(Pa{1},NVERT,latuface,longuface,'latitude (degree)','longitude (degree)','air pressure field (hpa)')

% Plot u
for i=1:Nt
plot3Dsliceresults(u{i},1,latuface,longuface,'latitude (degree)','longitude (degree)','u field (m/s)')
pause
clf
end 
% Plot v
plot3Dsliceresults(v{1},1,latvface,longvface,'latitude (degree)','longitude (degree)','v field (m/s)')






