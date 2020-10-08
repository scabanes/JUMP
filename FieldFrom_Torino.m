%########################################################################################################
%                                                                                           LOAD it -> 0:
%########################################################################################################
if(it==0)
%##########################################################################################
if(1==0) %###################################################################### DATA SET 1
% Here we load velocity fields:
[Vtheta_tot,b]=loadmtx([roots,Name,NameVt]);
[Vr_tot,b]=loadmtx([roots,Name,NameVr]);
% Here we map the loaded fields:
[returnOK] = Maps(Vtheta_tot,Vr_tot,Grid_Xp_cm,Grid_Yp_cm,GridR_2C,GridT_2C,Nti,Nri,Nrmin,Nrmax,Ntmin,Ntmax,theta,1)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            TRUNCATION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We truncated the domain size if necessary.
Grid_Xp_cm=Grid_Xp_cm(Ntmin:Ntmax,Nrmin:Nrmax);
Grid_Yp_cm=Grid_Yp_cm(Ntmin:Ntmax,Nrmin:Nrmax);
% Locally pseudo-cartesian grid see Read. 2015
GridR_2C=GridR_2C(Ntmin:Ntmax,Nrmin:Nrmax);
GridT_2C=GridT_2C(Ntmin:Ntmax,Nrmin:Nrmax);
% -----------------------------------------------------Truncation
theta = theta(1:length(Ntmin:Ntmax)); % Attenzione theta deve partire da zero..
r = r(Nrmin:Nrmax); % Attenzione theta deve partire da zero..
% Here we choose a pseudo cartesian grid 
x = Ro.*theta; 
y = r;
%######################################################################################
else %###################################################################### DATA SET 2
 file = [roots,Name,'/merged.civ2/fig_',num2str(numt),'.nc'];
% compter les images
numt = numt+gapt;
% Velocity fields of cartesian velocities
Vtheta = double(ncread(file,'U'));
Vr = double(ncread(file,'V'));
% regular cartesian coordinates
xx = double(ncread(file,'coord_x'));
yy = double(ncread(file,'coord_y'));
% Recreat the mesh grid
[XX,YY] = meshgrid(xx,yy);
[returnOK] = Maps_Cartesian(Vtheta,Vr,XX,YY,idxmin,idxmax,idymin,idymax)
end
%##########################################################################################
%##########################################################################################
end
%########################################################################################################
%                                                                                           LOAD it -> 1:
%########################################################################################################
if(it>0)
%##########################################################################################
%############################################################################### DATA SET 1
% Velocity fields of cartesian velocities
Vtheta = double(ncread(file,'U'));
Vr = double(ncread(file,'V'));
% regular cartesian coordinates
xx = double(ncread(file,'coord_x'));
yy = double(ncread(file,'coord_y'));
%
x = xx(idxmin:idxmax);
y = yy(idymin:idymax);
x=x-x(1); % we start at x=0
y=y-y(1); % we start at y=0
% beta
h = tan(11*pi/180)*[0.1:0.01:2.5]; % en m
dhdr = tan(11*pi/180); %
hbar=0.4; % en m
BETA = (2.*Omega/hbar) .* dhdr;
end