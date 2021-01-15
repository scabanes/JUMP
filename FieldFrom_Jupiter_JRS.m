%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mostly when derivative are involved, e.g energy and enstrophy fluxes.
% It is of prime interest to well define the sign of the velocity field on
% axis x and y as well as the direction of the derivatives along this axis
% by setting 'DevType_x' and 'DevType_y'. With Maps_Cartesian, we represent
% the velocity field as function of iy and ix indices. Look at the velocity
% field and give coherence with real field, then with the derivative.
%% Along the y axis, 
% the indices in the matrix are increasing downward on the vertical, if the 
% y velocity component is positive upward/downward then the derivative d/dy
%  must be descend/ascend along iy indices (it reverses axis and indices).
%% Along the x axis, 
% the indices in the matrix are increasing rightward on the horizontal, if  
% the x velocity component is positive rightward/leftward then the 
% derivative d/dy must be ascend/descend along iy indices (it reverses axis
% and idices).
%% => d/dx and d/dy must follow the direction of x and y axis given by the 
%%    positive Vx and Vy velocity component. 
%     1/ define the x and y axis: make your flow corresponds to its real
%     image, then positive velocity gives the direction of the axis.
%     2/ define the ascendent/descendent derivative along indices in order
%     to follow the x and y axis positive direction.
%########################################################################################################
%                                                                                           LOAD it -> 0:
%########################################################################################################
if(it==0)
% file to load    
file = [roots,Name,'/StatisticalData-dissip10000.nc'];
%---------------------------------------------------------------------LOAD:
% Velocity fields of cartesian velocities
Vtheta = double(ncread(file,'wm',[1,1,1,1],[361,720,1,1])); %wm(lat,long,z,time)
Vr = -double(ncread(file,'vm',[1,1,1,1],[361,720,1,1])); %vm(lat,long,z,time)
curl = double(ncread(file,'vort',[1,1,1,1],[361,720,1,1])); %vort(lat,long,z,time)
%
lat = double(ncread(file,'lat'))';
lon = double(ncread(file,'lon'))';
Colat =  (lat - 90).*sign(lat - 90).*(pi/180); % [0 pi]
%---------------------------------------------------------------------Grid:
yy = Colat.*R_Jup;
xx = lon.*(pi/180).*R_Jup;
% Recreat the mesh grid
[XX,YY] = meshgrid(xx,yy);
[returnOK] = Maps_Cartesian(Vtheta,Vr,XX,YY,idxmin,idxmax,idymin,idymax,DevType_x,DevType_y)
Ro=R_Jup;
end
%########################################################################################################
%                                                                                           LOAD it -> 1:
%########################################################################################################
if(it>0)
% Velocity fields of cartesian velocities
Vtheta = double(ncread(file,'wm',[1,1,1,it],[361,720,1,1])); %wm(lat,long,z,time)
Vr = -double(ncread(file,'vm',[1,1,1,it],[361,720,1,1])); %vm(lat,long,z,time)
curl = double(ncread(file,'vort',[1,1,1,it],[361,720,1,1])); %vort(lat,long,z,time)
% single time
% ----- cartesian coordinates
Ux = Vtheta;
Uy = Vr;
% curl = curl;
%
x = xx(idxmin:idxmax);
y = yy(idymin:idymax);
Nx = length(x);
Ny = length(y);
dx=mean(x(2:end)-x(1:end-1));
dy=mean(y(2:end)-y(1:end-1));
Lx = x(end)-x(1); % en m
Ly = y(end)-y(1); % en m
end