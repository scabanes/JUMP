function [returnOK] = Maps_Cartesian(Vx,Vy,XX,YY,idxmin,idxmax,idymin,idymax,DevType_x,DevType_y)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
addpath('/home/simon/Bureau/Esperimento-DICEA/JUMP/Altre-funzioni/cbrewer/')  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       General Mapped data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------Ux
fig=figure
imagesc(Vx(idymin:idymax,idxmin:idxmax));
colormap(cbrewer('div', 'RdBu', 50))
h1 = colorbar;
% caxis([-50 50])
set(get(h1,'label'),'string','sign(Vx) define x axis','FontSize',13); % sign(Vx) define x axis with respect to real flow
title('$V_x(iy,ix)$  -- truncated field','FontSize',18,'FontWeight','bold','Interpreter','Latex')
if(DevType_x == 'ascend_x') 
xlabel (['d/dx ascend on ix'],'FontSize',18,'FontWeight','bold','Interpreter','Latex')
elseif(DevType_x == 'descen_x') 
xlabel (['d/dx descend on ix'],'FontSize',18,'FontWeight','bold','Interpreter','Latex')
end
if(DevType_y == 'ascend_y') 
ylabel (['d/dy ascend on iy'],'FontSize',18,'FontWeight','bold','Interpreter','Latex')
elseif(DevType_y == 'descen_y') 
ylabel (['d/dy descend on iy'],'FontSize',18,'FontWeight','bold','Interpreter','Latex')
end
%-----------------------------------------Ux
fig=figure
imagesc(Vy(idymin:idymax,idxmin:idxmax));
colormap(cbrewer('div', 'RdBu', 50))
h2 = colorbar
% caxis([-50 50])
set(get(h2,'label'),'string','sign(Vy) define y axis','FontSize',13); % sign(Vx) define x axis with respect to real flow
title('$V_y(iy,ix)$  -- truncated field','FontSize',18,'FontWeight','bold','Interpreter','Latex')
if(DevType_x == 'ascend_x') 
xlabel (['d/dx ascend on ix'],'FontSize',18,'FontWeight','bold','Interpreter','Latex')
elseif(DevType_x == 'descen_x') 
xlabel (['d/dx descend on ix'],'FontSize',18,'FontWeight','bold','Interpreter','Latex')
end
if(DevType_y == 'ascend_y') 
ylabel (['d/dy ascend on iy'],'FontSize',18,'FontWeight','bold','Interpreter','Latex')
elseif(DevType_y == 'descen_y') 
ylabel (['d/dy descend on iy'],'FontSize',18,'FontWeight','bold','Interpreter','Latex')
end
%
% Ax = Vx(idymin:idymax,idxmin:idxmax);
% Ay = Vy(idymin:idymax,idxmin:idxmax);
% X = XX(idymin:idymax,idxmin:idxmax);
% Y = XX(idymin:idymax,idxmin:idxmax);
% figure
% quiver(flip(Ax(1:3:end,1:3:end)),Ay(1:3:end,1:3:end))
% set (gca,'Ydir','reverse')
% quiver(XX,YY,Vx,Vy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       General Mapped data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0==0)
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
%
returnOK = 'Ok';
end

