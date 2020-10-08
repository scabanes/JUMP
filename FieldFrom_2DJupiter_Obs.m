%########################################################################################################
%                                                                                           LOAD it -> 0:
%########################################################################################################
if(it==0)
% file to load    
file = [roots,Name,'/StatisticalData.nc'];
% Velocity fields of cartesian velocities
Vtheta_zt = double(ncread(file,'wm')); %wm(lat,long,z,time)
Vr_zt = double(ncread(file,'vm')); %vm(lat,long,z,time)
curl_zt = double(ncread(file,'vort')); %vort(lat,long,z,time)
% single time
Vtheta = Vtheta_zt(:,:,1,1);
Vr = Vr_zt(:,:,1,1);
curl = curl_zt(:,:,1,1);
%
SU=size(Vtheta);
% regular cartesian coordinates
xx=[1:(2.*pi.*R_Jup)./SU(2):2.*pi.*R_Jup];
yy=[1:(pi.*R_Jup)./SU(1):pi.*R_Jup];
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
Vtheta_zt = double(ncread(file,'wm')); %wm(lat,long,z,time)
Vr_zt = double(ncread(file,'vm')); %vm(lat,long,z,time)
curl_zt = double(ncread(file,'vort')); %vort(lat,long,z,time)
% single time
Vtheta = Vtheta_zt(:,:,1,it);
Vr = Vr_zt(:,:,1,it);
curl = curl_zt(:,:,1,it);
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