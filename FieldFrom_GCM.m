%########################################################################################################
%                                                                                           LOAD it -> 0:
%########################################################################################################
if(it==0)
% file to load    
file = [roots,Name,'/StatisticalData-SRS.nc'];
% Velocity fields of cartesian velocities
Vtheta = double(ncread(file,'wm',[1,1,1,1],[360,720,1,1])); %wm(lat,long,z,time)
Vr = double(ncread(file,'vm',[1,1,1,1],[360,720,1,1])); %vm(lat,long,z,time)
curl = double(ncread(file,'vort',[1,1,1,1],[360,720,1,1])); %vort(lat,long,z,time)
%
SU=size(Vtheta);
% regular cartesian coordinates
xx=[1:(2.*pi.*R)./SU(2):2.*pi.*R];
yy=[1:(pi.*R)./SU(1):pi.*R];
% Recreat the mesh grid
[XX,YY] = meshgrid(xx,yy);
[returnOK] = Maps_Cartesian(Vtheta,Vr,XX,YY,idxmin,idxmax,idymin,idymax)
Ro=R;
end
%########################################################################################################
%                                                                                           LOAD it -> 1:
%########################################################################################################
if(it>0)
% Velocity fields of cartesian velocities
Vtheta = double(ncread(file,'wm',[1,1,1,it],[360,720,1,1])); %wm(lat,long,z,time)
Vr = double(ncread(file,'vm',[1,1,1,it],[360,720,1,1])); %vm(lat,long,z,time)
curl = double(ncread(file,'vort',[1,1,1,it],[360,720,1,1])); %vort(lat,long,z,time)
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