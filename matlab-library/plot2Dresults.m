%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot2Dresults.m
% usage: plot 2D results from regional model
% only for cell centered 2D data
% Yun Zhang 05/07/2015
% @Stanford
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot2Dresults(input,x,y,xtype,ytype,titleinfo)
pcolor(x,y,input)
xlabel(xtype);
ylabel(ytype);
title(titleinfo);
colorbar;
colormap('Jet');
end
