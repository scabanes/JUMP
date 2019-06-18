%**********************************************************
%           PARAMETRI DA SCEGLIERE
%**********************************************************
%---------------------------------------------
%-------------- Nomenclature
%********************** Fields interpp on a polar grid 
NameVt = '/Vz29_60_360'; % Put the name of the azimuthal velocity file
NameVr = '/Vr29_60_360'; % Put the name of the radial velocity file
%********************** Fields interpp on a cartesian grid 
NameU = '/U2_cart';
NameV = '/V2_cart';
NameVort =  '/Vort2_cart';

%---------------------------------------------
%-------------- Fields dimensions
Nri=60; % number of points in r
Nti=360; % number of points in theta
Nrk=60; % number of modes in radius
nFrames = 2401; % number of Frames

%-----------------------------------------------------------
%-------------- Truncation parameters for Spectral analysis:
%-- on azimuth
%Ntmin=1; Ntmax=Nti;   % Keep full azimuthal domain
Ntmin=221; Ntmax=310;  % Narrow the domain in azimuth, !! for m interger
%-- on radius
Nrmin=1; Nrmax=Nri;    % Keep full radial domain 
%Nrmin=10; Nrmax=40;   % Narrow the domain in radius

%-----------------------------------------------------------
%-------------- Truncation parameters for PV monotonization:
PtDaTogliereR = [];% one has to put '[]' if no point to take out 
PtDaTogliereT = [1:Ntmin,Ntmax:360];% one has to put '[]' if no point to take out 

%---------------------------------------------
%-------------- Experimental device dimensions (cm)
g = 981;    % accelerazione di gravita
Lx = 68;    % dimensione orizzontale della vasca
Ly = 69;    % dimensione verticale della vasca
H0 = 4;     % altezza di fluido a riposo
Omega = 3;  % velocita angolare in rad/s

%---------------------------------------------
%-------------- Spectral quantities
Cz = 0.2;
Ck = 6;
beta = 0.53; % cm-1.s-1
epsilon = 2.1*10.^(-4); % cm-2.s-3

%---------------------------------------------
%-------------- Options
windowing = 0;  % if 0 no windowing is performed
ifsave = 1;     % if 0 no saving is performed
