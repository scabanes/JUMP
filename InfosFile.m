%**********************************************************
%           PARAMETRI DA SCEGLIERE
%**********************************************************
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
%-------------- Experimental device dimensions
R=29.7; % radius of the tank in cm

%---------------------------------------------
%-------------- Spectral quantities
Cz = 0.2;
Ck = 6;
beta = 0.53; % cm-1.s-1
epsilon = 2.1*10.^(-4); % cm-2.s-3