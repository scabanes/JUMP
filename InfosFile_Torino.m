%**********************************************************
%           PARAMETRI DA SCEGLIERE
%**********************************************************
% NOTE: Demander à Stefania de rajouter le nombre de frames 
% dans le titre du fichier des champs de vitesse, 'Vz29_60_360_nFrames'.
%---------------------------------------------
%-------------- Nomenclature
NameVt = '/Vz_540_180_1105'; % Put the name of the azimuthal velocity file
NameVr = '/Vr_540_180_1105'; % Put the name of the radial velocity file
%---------------------------------------------
%-------------- Frames
starti = 2;
starto = 1;
endi = 33122;
gapi = 1;
gapo = 2;
ftype = 'jpg';
step = 30;
%---------------------------------------------
%-------------- Fields dimensions
% La mapa Norm_nm permette di scegliere il numero di Nrk che condusce a un 
% erorre decente, i.e Norm_mn vicino a 1.
Nri=180; % number of points in r
Nti=540; % number of points in theta
Nrk=110; % number of modes in radius

%-----------------------------------------------------------
%-------------- Truncation parameters for Spectral analysis:
%->>> on azimuth_____________________
% Quando tagli un settore angolare, il settore deve essere 2*pi/alpha con  
% alpha un integer. Questo alpha define il modo azimutale m' il piu grande 
% come m'=alpha*m.
%
%Ntmin=1; Ntmax=Nti;   % Keep full azimuthal domain
Ntmin=103; Ntmax=156;  % alpha=10 ==> 2*pi/alpha  soit 54 points

%->>> on radius______________________
% Possiamo tagliare in ma se togliamo in primi punti tra r=[0:rT]<R
% ma i raggi non possono essere r=[rT:R] ma invece devono essere ridefiniti
% come r=[0:R-rT] e R'=R-rT. Cosi rispettiamo le proprieta delle funzioni
% di Bessel, i.e Norm_mn +-= 1. Facendo questo l integrale rispetta:
% Int_0^R' r*u(m,r)*Jm(alpha_mn*r/R') dr con r=[0:R-rT] <=> r=[0:R'] ; 
%
%Nrmin=1; Nrmax=Nri;      % Keep full radial domain 
Nrmin=48; Nrmax=Nri-10;   % Narrow the domain in radius

%-----------------------------------------------------------
%-------------- Truncation parameters for PV monotonization:
PtDaTogliereR = [1, 178, 179,180];%[];% one has to put '[]' if no point to take out 
PtDaTogliereT = [];%[1:Ntmin,Ntmax:360];% one has to put '[]' if no point to take out 

%---------------------------------------------
%-------------- Experimental device dimensions (cm)
g = 981;    % accelerazione di gravita
H0 = 56;     % altezza di fluido a riposo
Lx = 500;    % dimensione orizzontale della vasca
Ly = 500;    % dimensione verticale della vasca
Omega = 0.16;  % velocita angolare in rad/s

%---------------------------------------------
%-------------- Spectral quantities
Cz = 0.5;
Ck = 6;
beta = 0.11/100; % cm-1.s-1
epsilon = 1.*10.^(-5); % cm-2.s-3

%---------------------------------------------
%-------------- Options
windowing = 1;  % if 0 no windowing is performed
ifsave = 1;     % if 0 no saving is performed



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncation parameters
% Ntmin=1; Ntmax=Nti;%Nt-50;%180; %limits in theta and r, if you need to narrow the domain
% Ntmin=79; Ntmax=168; % alpha=6 ==> 2*pi/alpha soit 90 points
% Ntmin=100; Ntmax=159; % alpha=9 ==> 2*pi/alpha  soit 60 points
% Ntmin=103; Ntmax=156; % alpha=10 ==> 2*pi/alpha  soit 54 points
