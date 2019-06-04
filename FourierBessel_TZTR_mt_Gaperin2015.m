clear all
close all
%**********************************************************
%           PARAMETRI DA SCEGLIERE
%**********************************************************
%roots = 'D:/simon/esperimento30/';
roots = 'D:\Gabriella\ANALISI\LUGLIO 2015\ANALISI DAI FES\Jesse\1 Arco  90°\analisi da fes\Ncerchi60_Nraggi120\esp 30/'
% Resoluzione dei dati
Nri=60; % number of points in r
Nti=120; % number of points in theta
Nrk=60; % number of modes in radius
%
Create_Grid_pol_OK
% Truncation parameters
Ntmin=1; Ntmax=Nti;%Nt-50;%180; %limits in theta and r, if you need to narrow the domain
%Ntmin=221; Ntmax=310;
Ntmin=75; Ntmax=104;
Nrmin=1; Nrmax=Nri;
% Dimensioni
R=29.7; % radius of the tank in cm
scaleR=R/Nri;
r=(1:Nri)*scaleR;
theta = [0:(2*pi)/Nti:2*pi-(2*pi)/Nti];
%**********************************************************
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                   LOAD:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('D:\simon\script\besselzeros2_C.mat'); 
% -----------------------------------------------------Campo di velocita
g=2;
% [Vtheta_tot,a] = loadmtx(['..\Vz' num2str(g)]);
%[Vtheta_tot,a] = loadmtx([roots,'Vz30_360']);
[Vtheta_tot,a] = loadmtx([roots,'Vz2']);

%
%[Vr_tot,a] = loadmtx([roots,'Vr30_360']);
% [Vr_tot,a] = loadmtx(['..\Vr' num2str(g)]);
[Vr_tot,a] = loadmtx([roots,'Vr2']);

% -----------------------------------------------------Truncation
theta = theta(1:length(Ntmin:Ntmax)); % Attenzione theta deve partire da zero..
% !!!!!! Attention encore un doute
%r = r(1:length(Nrmin:Nrmax)); % Attenzione theta deve partire da zero..
%R=r(end);
%
% r = r(Nrmin:Nrmax); % Attenzione theta deve partire da zero..
% !!!!!! -------------------------
Nt = length(theta);
Nr = length(r);
%

Tmax=3;%a(2);
EZn = zeros(Nrk,Tmax);
for it=1:Tmax
    
Vtheta = reshape(Vtheta_tot(:,it),Nti,Nri);
Vr = reshape(Vr_tot(:,it),Nti,Nri);
Vtheta=Vtheta(Ntmin:Ntmax,Nrmin:Nrmax);
Vr=Vr(Ntmin:Ntmax,Nrmin:Nrmax);

if(0==1)
Grid_Xp_cm=Grid_Xp_cm(Ntmin:Ntmax,Nrmin:Nrmax);
Grid_Yp_cm=Grid_Yp_cm(Ntmin:Ntmax,Nrmin:Nrmax);

figure; title('section we are working on')
contourf(Grid_Xp_cm,Grid_Yp_cm,Vtheta)
end
%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   DECOMPOSITION IN SPECTRAL DOMAIN
%  ###########################################################################################################################
%  ###########################################################################################################################

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       FOURIER DECOMPOSITION IN AZIMUTH:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% la fft deve essere normalisata per il numero di punti in azimuth che corrisponde
% al numero dei modi Nm che escono della fft. Tra questi Nm modi ci sono M modi
% zonali m che corrisponde al numero dei punti unici simetrici per M ne l vettore _fft.
%-----------Vtheta
Vtheta_fft=fft(Vtheta,[],1)/size(Vtheta,1); % --> Vtheta_fft(m,r)
%-----------Vr
Vr_fft=fft(Vr,[],1)/size(Vr,1); % --> Vtheta_fft(m,r)
%-----------modes
Nm = size(Vtheta_fft,1);
M=ceil((size(Vtheta_fft,1)+1)/2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         BESSEL DECOMPOSITION IN RADIUS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Troviamo qui i coefficienti "Cmn" usando la formula del libro Arfken et al. 
% Mathematical Methods for Physicists, chapter 14 Bessel function, expression
% (14.48). Questa formula verifica la normalisazione expression (14.46) che
% cacoliamo come "Norm_mn" e deve essere uguale a 1 (al meno essere vicino).
Cmn_U = zeros(Nm,Nrk); % -->  Cmn_U(m,n)
Cmn_V = zeros(Nm,Nrk); % -->  Cmn_V(m,n)
vm=zeros(Nm,1); % --> vm(m)
Norm_mn = zeros(Nm,Nrk);
% for positive frequences m
for m=1:M % indice matriciel, m-1 pour le mode
    for nr=1:Nrk
     %-----------Vtheta
     integrand=r.*besselj(m-1,J_root(m,nr).*r/R).*(Vtheta_fft(m,:));
     Cmn_U(m,nr)=2./((R*besselj(m,J_root(m,nr))).^2).*trapz(r,integrand);
     %-----------Vr
     integrand=r.*besselj(m-1,J_root(m,nr).*r/R).*(Vr_fft(m,:));
     Cmn_V(m,nr)=2./((R*besselj(m,J_root(m,nr))).^2).*trapz(r,integrand);
     %-----------Norm_mn
     integrand_N=r.*besselj(m-1,J_root(m,nr).*r/R).*besselj(m-1,J_root(m,nr).*r/R);
     Norm_mn(m,nr)=2./((R*besselj(m,J_root(m,nr))).^2).*trapz(r,integrand_N);
    end
    vm(m)=m-1; %vecteur des modes m
end
% for negative frequences m
mm=0;
for m=Nm:-1:M+1 % indice matriciel
   mm=mm+1;     % mm+1 indice matriciel pour les zeros de bessel, -mm pour le mode
   for nr=1:Nrk
     %-----------Vtheta
     integrand=r.*besselj(-mm,J_root(mm+1,nr).*r/R).*(Vtheta_fft(m,:));
     Cmn_U(m,nr)=2./((R*besselj(-mm-1,J_root(mm+1,nr))).^2).*trapz(r,integrand);
     %-----------Vr
     integrand=r.*besselj(-mm,J_root(mm+1,nr).*r/R).*(Vr_fft(m,:));
     Cmn_V(m,nr)=2./((R*besselj(-mm-1,J_root(mm+1,nr))).^2).*trapz(r,integrand);
     %-----------Norm_mn
     integrand_N=r.*besselj(-mm,J_root(mm+1,nr).*r/R).*besselj(-mm,J_root(mm+1,nr).*r/R);
     Norm_mn(m,nr)=2./((R*besselj(-mm-1,J_root(mm+1,nr))).^2).*trapz(r,integrand_N);
   end
   vm(m)=-mm;  %vecteur des modes m
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                SPECTRA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ab2_temp = Cmn_U.*conj(Cmn_U) + Cmn_V.*conj(Cmn_V);
Emn = zeros(Nm,Nrk);
% for positive frequences m
for m=1:M % indice matriciel, m-1 pour le mode
Emn(m,:) = 0.5.*ab2_temp(m,:).*besselj(m,J_root(m,1:Nrk)).^2.;
end
% for negative frequences m
mm=0;
for m=Nm:-1:M+1 % indice matriciel
   mm=mm+1;     % mm+1 indice matriciel pour les zeros de bessel, -mm pour le mode
Emn(m,:) = 0.5.*ab2_temp(m,:).*besselj(-mm-1,J_root(mm+1,1:Nrk)).^2.;
end

EZn(:,it) = Emn(1,:);
ERn(:,it) = sum(Emn(2:end,:));
end
%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   PROJECTION IN SPACIAL DOMAIN
%  ###########################################################################################################################
%  ###########################################################################################################################

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                   PROJECTION ON BESSEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vtheta_fft_proj=zeros(Nm,Nr);
% for positive frequences m
for m=1:M%Nt2
for ir=1:Nr
Vtheta_fft_proj(m,ir) = sum( Cmn_U(m,:) .* besselj(m-1,J_root(m,1:Nrk).*r(ir)/R));
end
end
% for negative frequences m
mm=0;
for m=Nm:-1:M+1
   mm=mm+1;
for ir=1:Nr
Vtheta_fft_proj(m,ir) = sum( Cmn_U(m,:) .* besselj(-mm,J_root(mm+1,1:Nrk).*r(ir)/R));
end 
end
% -------------------------------- plots
figure; hold on
plot(real(Vtheta_fft(20,:)))
plot(real(Vtheta_fft_proj(20,:)),'*')
title('Reconstruct I')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                  PROJECTION ON FOURIER:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% costruiamo il vettore m che corrisponde anche al vettore vm. Era solo per
% veriificare.
clear m
m = zeros(Nm,1);
m(1:M) = [1:M]-1;
m(M+1:Nm)=-(M-2:-1:1);
% Se tagliamo in theta i modi diventano m = alpha*m
% alpha e la frazione del settore angolare.
alpha=round((2*pi)/(theta(end)-theta(1)))
disp(['Attenzione m = ',num2str(alpha),'*m'])
%
Vtheta_proj= zeros(size(Vtheta));
Vtheta_proj_step2= zeros(size(Vtheta));
% loop sui raggi e theta per rimappare
for ir=1:Nr
    for itheta=1:Nt
       % Vtheta_proj(itheta,ir)= sum(Vtheta_fft(:,ir) .* exp((-2.*pi*i)/M.*-(m-1).*(itheta-1))); 
        Vtheta_proj(itheta,ir)= sum(Vtheta_fft(:,ir) .* exp((1i).*(alpha.*m).*theta(itheta))); 
        Vtheta_proj_step2(itheta,ir)= sum(Vtheta_fft_proj(:,ir) .* exp((1i).*(alpha.*m).*theta(itheta))); 
    end
end
% -------------------------------- plots maps
if(0==1)
figure; hold on
contourf(Vtheta)
figure; hold on
contourf(real(Vtheta_proj))
figure; hold on
contourf(real(Vtheta_proj_step2))
end
% -------------------------------- plots profils
figure
plot(real(Vtheta(:,25)))
hold on
plot(real(Vtheta_proj(:,25)),'*')
plot(real(Vtheta_proj_step2(:,25)),'-*')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                SPECTRA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -------------------------------- Theoretical spectra
Cz = 0.5;
Ck = 6;
beta = 0.53; % cm-1.s-1
epsilon = 11*10.^(-4); % cm-2.s-3
EZ_Theo = (Cz*beta.^2.).*(1./R).*(J_root(1,1:Nrk)./R).^(-5.);
ER_Theo = Ck.*(epsilon.^(2./3.))*(1./R).*(J_root(2,1:Nrk)./R).^(-5./3.);
n = [1:Nrk]; % non-dimensional
% -------------------------------- plots spectra
figure
loglog(n,mean(EZn,2),'r')
hold on
loglog(n,mean(ERn,2),'k')
loglog(n,EZ_Theo,'--r')
loglog(n,ER_Theo,'--k')
ylim([10.^(-7.) 0.1])
xlim([n(1) n(end)])


Ispectral = sum(mean(EZn,2)) + sum(mean(ERn,2));
disp('####################################')
disp('Energy = '), disp(num2str(Ispectral))
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                      SURFACE INTEGRALS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = R*R*2*pi/(2.*alpha);%pi*R*R; % non serve qui per gli spectri, la superficie scompare con l integrale su theta
% ForInt = zeros(size(Vtheta));
for itheta=1:Nt
    IAnalytique_rT = r.* 0.5.* (abs(Vtheta(itheta,:)).^2 + abs(Vr(itheta,:)).^2);
    IAnalytique_T(itheta) = trapz(r,IAnalytique_rT);
end
IAnalytique = trapz(theta,IAnalytique_T)./A;

% I = trapz(y,trapz(x,F,2))

    
disp('####################################')
disp('Energy = '), disp(num2str(IAnalytique))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for itheta=1:Nt
    IA_rT = r.* 0.5.* abs(Vtheta_proj_step2(itheta,:)).^2;
    IA_T(itheta) = trapz(r,IA_rT);
end
IA = trapz(theta,IA_T)./A;

% I = trapz(y,trapz(x,F,2))
   
disp('####################################')
disp('Energy = '), disp(num2str(IA))

