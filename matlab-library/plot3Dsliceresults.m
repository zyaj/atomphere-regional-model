%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot3Dsliceresults.m
% usage: plot 3D results for a specific layer from regional model
% Yun Zhang 05/07/2015
% @Stanford
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot3Dsliceresults(input,k,x,y,xtype,ytype,titleinfo)
[a,b]=size(x);
plotdata=reshape(input(:,k),b,a);
plotdata=plotdata';
pcolor(x,y,plotdata)
xlabel(xtype);
ylabel(ytype);
title(titleinfo);
colorbar;
colormap('Jet');
end