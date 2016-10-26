%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotvelocityquiver.m
% usage: plot velocity quiver for specific layer
% Yun Zhang 05/07/2015
% @Stanford
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotvelocityquiver(u,v,k,x,y,xtype,ytype,titleinfo)
[a,b]=size(x);
uu=reshape(u(:,k),b+1,a);
vv=reshape(v(:,k),b,a+1);
vv=vv';
uu=uu';
uc=(uu(:,1:b)+uu(:,2:b+1))/2;
vc=(vv(1:a,:)+vv(2:a+1,:))/2;
quiver(x,y,uc,vc,1.4)
xlabel(xtype);
ylabel(ytype);
title(titleinfo);
end