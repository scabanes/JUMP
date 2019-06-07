%**********************************************************
%           PARAMETRI DA SCEGLIERE
%**********************************************************
% NOTE: Demander à Stefania de rajouter le nombre de frames 
% dans le titre du fichier des champs de vitesse, 'Vz29_60_360_nFrames'.
%---------------------------------------------
%-------------- Nomenclature
NameVt = 'Vz29_60_360'; % Put the name of the azimuthal velocity file
NameVr = 'Vr29_60_360'; % Put the name of the radial velocity file

%---------------------------------------------
%-------------- Fields dimensions
Nri=60; % number of points in r
Nti=360; % number of points in theta
Nrk=60; % number of modes in radius

%---------------------------------------------
%-------------- Truncation parameters
%-- on azimuth
%Ntmin=1; Ntmax=Nti;   % Keep full azimuthal domain
Ntmin=221; Ntmax=310;  % Narrow the domain in azimuth, !! for m interger
%-- on radius
Nrmin=1; Nrmax=Nri;    % Keep full radial domain 
%Nrmin=10; Nrmax=40;   % Narrow the domain in radius

%---------------------------------------------
%-------------- Frames wanted
% Tmax=nan; % if nan inputs the default value will be the number of frames 
        % contain in the input velocity file.
Tmax=1000;
Itime=1;

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
