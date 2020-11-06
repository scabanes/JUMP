function [returnOK] = Maps_Cartesian(Vx,Vy,XX,YY,idxmin,idxmax,idymin,idymax)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
addpath('/home/simon/Bureau/Esperimento-DICEA/JUMP/Altre-funzioni/cbrewer/')  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       General Mapped data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1==0)
%-----------------------------------------Ux
fig=figure
pcolor(XX,YY,Vx); shading interp
colormap(cbrewer('div', 'RdBu', 50))
colorbar
mycmap = get(fig,'Colormap')
set(fig,'Colormap',flipud(mycmap))
% caxis([-3 3]);
hold on
scatter(XX(idymin,idxmin),YY(idymin,idxmin), [], 'Xr')
scatter(XX(idymin,idxmax),YY(idymin,idxmax), [], 'Xr')
scatter(XX(idymax,idxmin),YY(idymax,idxmin), [], 'Xr')
scatter(XX(idymax,idxmax),YY(idymax,idxmax), [], 'Xr')
title('Full field and truncation: Vx')
%-----------------------------------------Uy
fig=figure
pcolor(XX,YY,Vy); shading interp
colormap(cbrewer('div', 'RdBu', 50))
colorbar
mycmap = get(fig,'Colormap')
set(fig,'Colormap',flipud(mycmap))
% caxis([-3 3]);
hold on
scatter(XX(idymin,idxmin),YY(idymin,idxmin), [], 'Xr')
scatter(XX(idymin,idxmax),YY(idymin,idxmax), [], 'Xr')
scatter(XX(idymax,idxmin),YY(idymax,idxmin), [], 'Xr')
scatter(XX(idymax,idxmax),YY(idymax,idxmax), [], 'Xr')
title('Full field and truncation: Vy')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       General Mapped data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------Ux
fig=figure
imagesc(Vx(idymin:idymax,idxmin:idxmax));
colormap(cbrewer('div', 'RdBu', 50))
colorbar
title('Vx truncated')
%-----------------------------------------Ux
fig=figure
imagesc(Vy(idymin:idymax,idxmin:idxmax));
colormap(cbrewer('div', 'RdBu', 50))
colorbar
title('Vy truncated')
%
Ax = Vx(idymin:idymax,idxmin:idxmax);
Ay = Vy(idymin:idymax,idxmin:idxmax);
X = XX(idymin:idymax,idxmin:idxmax);
Y = XX(idymin:idymax,idxmin:idxmax);
figure
quiver(Ax(1:3:end,1:3:end),Ay(1:3:end,1:3:end))
set (gca,'Ydir','reverse')
% quiver(XX,YY,Vx,Vy)
%
returnOK = 'Ok';
end

