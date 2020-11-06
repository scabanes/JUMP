%########################################################################################################
%                                                                                           LOAD it -> 0:
%########################################################################################################
if(it==0)
% file to load    
file = [roots,Name,'/StatisticalData.nc'];
%---------------------------------------------------------------------LOAD:
% Velocity fields of cartesian velocities
Vtheta = double(ncread(file,'wm',[1,1,1,1],[361,720,1,1])); %wm(lat,long,z,time)
Vr = double(ncread(file,'vm',[1,1,1,1],[361,720,1,1])); %vm(lat,long,z,time)
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
[returnOK] = Maps_Cartesian(Vtheta,Vr,XX,YY,idxmin,idxmax,idymin,idymax)
Ro=R_Jup;
end
%########################################################################################################
%                                                                                           LOAD it -> 1:
%########################################################################################################
if(it>0)
% Velocity fields of cartesian velocities
Vtheta = double(ncread(file,'wm',[1,1,1,it],[361,720,1,1])); %wm(lat,long,z,time)
Vr = double(ncread(file,'vm',[1,1,1,it],[361,720,1,1])); %vm(lat,long,z,time)
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