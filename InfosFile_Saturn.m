%**********************************************************
%           PARAMETRI DA SCEGLIERE
%**********************************************************
%---------------------------------------------
%-------------- Nomenclature
%********************** Fields interpp on a polar grid 
NameVt = '/UPN'; % Put the name of the azimuthal velocity file
NameVr = '/VPN'; % Put the name of the radial velocity file
%********************** Fields interpp on a cartesian grid 
%NameU = '/Ux_ee_PN_cart'; % in the x direction - eddy-eddy velocity (-zonal mean of the zonal velocity)
%NameV = '/Uy_ee_PN_cart'; % in the y direction - eddy-eddy velocity (-zonal mean of the zonal velocity)
NameU = '/UxPN_cart'; % in the x direction - total velocity
NameV = '/UyPN_cart'; % in the y direction - total velocity
NameVort =  '/vortPN_cart';

%---------------------------------------------
load(['Data2RegularGrid_infos.mat'])
%-------------- Fields dimensions
%Nri=133; % number of points in r
%Nti=200; % number of points in theta
Nrk=60; % number of modes in radius
nFrames = 1; % number of Frames

%-----------------------------------------------------------
%-------------- Truncation parameters for Spectral analysis:
%-- on azimuth
Ntmin=1; Ntmax=Nti;   % Keep full azimuthal domain
%Ntmin=1; Ntmax=14;  % Narrow the domain in azimuth, !! for m interger
%-- on radius
Nrmin=1; Nrmax=Nri;    % Keep full radial domain 
%Nrmin=10; Nrmax=40;   % Narrow the domain in radius

%-----------------------------------------------------------
%-------------- Truncation parameters for PV monotonization:
% qui si deve mettere i punti da togliere
PtDaTogliereR = [];% one has to put '[]' if no point to take out 
PtDaTogliereT = [];%[1:Ntmin,Ntmax:360];% one has to put '[]' if no point to take out 

%---------------------------------------------
%-------------- Experimental device dimensions (cm)
% % g = 981;    % accelerazione di gravita
% % Lx = 68;    % dimensione orizzontale della vasca
% % Ly = 69;    % dimensione verticale della vasca
% % H0 = 4;     % altezza di fluido a riposo
Omega = 0.000165121;  % velocita angolare in rad/s
R_sat = 58232000; % en m

%---------------------------------------------
%-------------- Spectral quantities
Cz = 0.2;
Ck = 6;
% beta = 10.^(-12); % cm-1.s-1
beta = ((2*Omega)/R_sat)*cos(2*pi/5);%% latitude 72° -- cos(pi/3); latitude 60°
epsilon = 4.*10.^(-5); % m2.s-3
%epsilon = 2.*10.^(-6); % m2.s-3
nhat=15;

%---------------------------------------------
%-------------- Options
windowing = 0;  % if 0 no windowing is performed
ifsave = 1;     % if 0 no saving is performed


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					      PTS: parameter to set for EFluxes_Spectral.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is the number of l in pixels on which I calculate the fluxes. It is 
% theoretically limited by Nx and Ny the size of the velocity field in 
% pixels. But it appears that with Nl = Nx & Ny, if the velocity field is
% not periodic the filtering does not manage to filter out every large scales 
% structures. If you impose Nl >> Nx then you finally get to a flux = 0 at 
% largest scales. This might be improved by using a different filter than a
% currently used gaussian, e.g une fonction porte. This can be changed in
% func_Spatial_GaussianFilter.m and/or func_Spectral_GaussianFilter.m
Nl = 120;%128;
%--------------------------------------------
%		Spatial average area
%--------------------------------------------
% This is truncating the size of the domain we compute the averaged fluxes on. 
% It is choosen to truncate the edge of the domain where fluxes are artificially
% strong due to the non periodicity of the velocity maps.
%idxtc=60; % This is the truncation parameter on x
%idytc=100; % This is the truncation parameter on y
%--------------
if(0==0)
xmin=6%;1;
xmax=114;%120;
ymin=6;%1;
ymax=114;%120;
radius = 0; % if(0) no circle, if(1..) croop of a centred circle with a radius in pixels.
else
xmin=30;
xmax=100;
ymin=1;
ymax=60;	
end
%--------------------------------------------
% In spectral filtering it results that l=1 is an uresolved scale that
% corresponds to the lunquist frequency and is theoretically out of
% resolution, consequently the resulting filtered flow is meaning less at l=1. 
% Spatial filtering however manage with it by filtering nothing.
vl=[1:1:Nl];%([0.5:1:Nl]);
% Here we can c
FSpatial=0;% set to 1 to have spatial filter and 0 to have spectral filter.
gaussian = 2; % If gaussian = 1: then G~exp(-6r^2/l^2) if 2: G~exp(-9r^2/(2l^2))
%--------------------------------------------
%  Spatial cartesian area to compute fluxes
%--------------------------------------------
% In EFluxes_Spectral.m I work on a squared 
% cartesian grid:
%---------------------------------------------
% Cartesien Nx=Ny=
if(0==0)
 idxmin = 1; %33
 idxmax = 120; %96
 idymin = 1;
 idymax = 120; %64
else
 idxmin = 1;
 idxmax = 128;
 idymin = 1;
 idymax = 128;
end
