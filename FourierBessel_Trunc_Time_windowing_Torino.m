%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code leads to azimuthal fourier decomposition and radial bessel
% decomposition of 2D field defined on a polar grid of Nri radial points and
% Nti azimuthal points. One obtains the complexe coefficients of the zonal 
% velocity Vtheta_tot denoted Cmn_U and radial veolcity Vr_tot denoted
% Cmn_V amd the axisymmetric and residual energy as a function of radial
% modes n and time.
% _________________________________________________________________________
% Truncation options:
% 1) Azimuthal truncation of the field, is obtained by giving truncated
% indices Ntmin and Ntmax. !!! Attention !!! These indices define also the
% new array of theta on which the fourier decomposition is performed. Theta
% goes as theta = [0:2*pi/alpha] and alpha has to be an integer which define
% the largest wave length in azimuth as exp(i*alpha*m*theta) with the
% azimuthal modes m=[0:M-1] for positive frequences and m=[1:M-2] for
% negative frequences. !!! Attention !!! It is then necessary to choose 
% Ntmin and Ntmax leading to an angular sector going from 0 to 2*pi/alpha.
% Note that the number of zonal modes M is defined by the number of
% theta points Nti once truncated, Nt, such that M=ceil(Nt+1)/2 ;
% 2) Radial truncation of the field, is obtained by giving truncated
% indices Nrmin and Nrmax. Such a truncation will reduce the number of radial 
% frequencies Nrk which define the radial bessel modes n = [1:Nkr] and the 
% radius array r=[0:R_trunc]. The total radius R_trunc is total radius after
% radial truncation. For the bessel function to be correctly normalised, r 
% as to be defined from 0 to R.
% The trucation in radius is recomanded when a region is not covered by
% data. In such a region, zero values will reduce the total energy in the
% spectra and will strongly affect the inversion fit to obtain
% Bessel-Fourier coefficients. Small zero regions are however tolerate.
% _________________________________________________________________________
% Normalisation to set Nrk:
% To have well defined radial modes we constrain our choice for 'Nrk' by
% looking at the normalization matrice 'Norm_mn' defined in  "Arfken et al. 
% Mathematical Methods for Physicists", expression (14.46). Matrice
% 'Norm_mn' has to be approximatively unity. If such is not the case one has
% to reduce Nrk.
% _________________________________________________________________________
% Windowing option:
% windowing = 1 is active & inactive with windowing = 0.
% The windowing in radial and azimuthal direction is performed via function 
% WR and WA that can use Hanning or Tukey function with parameter alpha
% within the windowing fuction that fix the extension of the smoothing. It
% sometimes help to improve the fit when inverting for BF coefficients.
% Note that when a 2*pi data set is avalaible, the azimuthal signal is 
% automatically periodic and azimuthal windowing has to be removed. This
% step is not mandatory...
% _________________________________________________________________________

% INPUTS: 
% # Set PDS parameters...
% # 'Vtheta_tot' and 'Vr_tot' the 2D velocity fields of azimuthal and
%    radial velocity of size (Nti*Nri,time).
% # Nri and Nti radial and azumuthal point numbers before truncation.
% # Ntmin and Ntmax are truncated azimuthal indices.
% # Nrmin and Nrmax are truncated radial indices.

% OUTPUTS:
% # Nr and Nt radial and azumuthal point numbers after truncation.
% # R_trunc, the maximum radius after truncation.
% # EZn, energy in the axisymmetric mode as a function of radial modes n
%   and time.
% # ERn, energy in the residual as a function of radial modes n and time.
% # Norm_mn, the normalisation matrice that help to choose Nrk.
% # error_fit, is the RMS of the difference between the velocity field
%   input and the one obtains from Bessel Fourier projection from
%   coefficients Cmn_U and Cmn_V.
% _________________________________________________________________________

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%**********************************************************
%          PDS: PARAMETRI DA SCEGLIERE
%**********************************************************
cexp = 6;
expt = {'5';'6';'7';'8';'9';'3'}
roots = '../Dati-Torino/';
output = 'outputs'
windowing = 1; % if 0 no windowing is performed
ifsave = 1; % if 0 no saving is performed
%**********************************************************
% ##################################################################################################################################################
%                                                                                                                                 DATI INFORMAZIONI:
% ##################################################################################################################################################
EXPT_infos
load([roots,names{cexp,2},'/',output,'/infos.mat'])
nFrames = ((endi{cexp,2}-starti{cexp,2})/step + 1);
Tmax=500;
load('D:\simon\script\besselzeros2_C.mat'); 
% Resoluzione dei dati
Create_Grid_pol_Torino_cm
Nri=Ncerchi; % number of points in r
Nti=Nraggi; % number of points in theta
% La mapa Norm_nm permette di
% scegliere il numero di Nrk che condusce a un erorre decente, i.e Norm_mn 
% vicino a 1.
Nrk=125%137; % number of modes in radius
Nrk=110%137; % number of modes in radius
% Truncation parameters
%%%%%%%%%%%%%%
% Quando tagli un settore angolare, il settore deve essere 2*pi/alpha con alpha un 
% integer. Questo alpha define il modo azimutale m' il piu grande come
% m'=alpha*m.
Ntmin=1; Ntmax=Nti;%Nt-50;%180; %limits in theta and r, if you need to narrow the domain
Ntmin=79; Ntmax=168; % alpha=6 ==> 2*pi/alpha soit 90 points
%Ntmin=87; Ntmax=164; % alpha=7 ==> 2*pi/alpha soit 77.1 points
Ntmin=100; Ntmax=159; % alpha=9 ==> 2*pi/alpha  soit 60 points
Ntmin=103; Ntmax=156; % alpha=10 ==> 2*pi/alpha  soit 54 points
%%%%%%%%%%%%%%
% POssiamo tagliare in ma se togliamo in primi punti tra r=[0:rT]<R
% ma i raggi non possono essere r=[rT:R] ma invece devono essere ridefiniti
% come r=[0:R-rT] e R'=R-rT. Cosi rispettiamo le proprietà delle funzioni
% di Bessel, i.e Norm_mn +-= 1. Facendo questo l integrale rispetta:
% Int_0^R' r*u(m,r)*Jm(alpha_mn*r/R') dr con r=[0:R-rT] <=> r=[0:R'] ; 
Nrmin=1; Nrmax=Nri;
Nrmin=42; Nrmax=Nri;
Nrmin=48; Nrmax=Nri-10;
%%%%%%%%%%%%%%
% Dimensioni
R=MaxRadius;%29.7; % radius of the tank in cm
scaleR=R/Nri;
r=(1:Nri)*scaleR;
theta = [0:(2*pi)/Nti:2*pi-(2*pi)/Nti];
% Grid to plot maps
% nFrames=2150;
% ##################################################################################################################################################
%                                                                                                                                              LOAD:
% ##################################################################################################################################################
% Carichiamo input--------------------------
name_mtx1 = [roots,names{cexp,2},'/',output,'/Vz_',num2str(Nraggi),'_',num2str(Ncerchi),'_',num2str(nFrames)];
[Vtheta_tot,b]=loadmtx(name_mtx1);
name_mtx2 = [roots,names{cexp,2},'/',output,'/Vr_',num2str(Nraggi),'_',num2str(Ncerchi),'_',num2str(nFrames)];
[Vr_tot,b]=loadmtx(name_mtx2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                             TRUNCATION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Grid_Xp_cm=Grid_Xp_cm(Ntmin:Ntmax,Nrmin:Nrmax);
Grid_Yp_cm=Grid_Yp_cm(Ntmin:Ntmax,Nrmin:Nrmax);
% -----------------------------------------------------Truncation
theta = theta(1:length(Ntmin:Ntmax)); % Attenzione theta deve partire da zero..
% !!!!!! Attention encore un doute
r = r(1:length(Nrmin:Nrmax)); % Attenzione theta deve partire da zero..
R=r(end);
% r = r(Nrmin:Nrmax); % Attenzione theta deve partire da zero..
% !!!!!! -------------------------
Nt = length(theta)
Nr = length(r)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     WINDOWING MATRICES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tukey or Hanning windowing are possible by respectively
% Windowing.WindowingTukey and Windowing.WindowingHanning. One has to fixe
% a alpha velue within the windowing function to reduce or increase the
% windowing extension to smooth to zero at both side.
WA = Windowing.WindowingTukey(Grid_Xp_cm,1); % Tukey function on the Azimuthal direction.
WR = Windowing.WindowingTukey(Grid_Xp_cm,2); % Tukey function on the Radial direction.
%%%%%%%%%%%%%%%%%%%%%%%ZEROS MATRICES
EZn = zeros(Nrk,Tmax);
ERn = zeros(Nrk,Tmax);
%  ###########################################################################################################################
%  ##################################################### Time Loop ###########################################################
%  ###########################################################################################################################
for it=1:Tmax
% it=Tmax;
     disp(['step = ',num2str(it),' on ',num2str(Tmax)])
%----------------- Reshape
Vtheta = reshape(Vtheta_tot(:,it),Nti,Nri);
Vr = reshape(Vr_tot(:,it),Nti,Nri);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(0==0)
% figure; title('section we are working on')
% contourf(Grid_Xp_cm,Grid_Yp_cm,Vtheta*100,15)%,'LineStyle','none')
% colorbar
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------- Truncation & mise en cm
Vtheta=Vtheta(Ntmin:Ntmax,Nrmin:Nrmax)*100;
Vr=Vr(Ntmin:Ntmax,Nrmin:Nrmax)*100;
%----------------- Isnan = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(0==0)
% figure; title('section we are working on')
% contourf(Grid_Xp_cm,Grid_Yp_cm,Vtheta,15)%,'LineStyle','none')
% colorbar
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vtheta(isnan(Vtheta))=0.;
Vr(isnan(Vr))=0.;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              WINDOWING:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(windowing==1) 
% ------------------------Azimuthal windowing
Vtheta = Vtheta.*WA;
Vr = Vr.*WA;
% ------------------------Radial windowing
Vtheta = Vtheta.*WR;
Vr = Vr.*WR;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(0==0)
% figure; title('section we are working on')
% contourf(Grid_Xp_cm,Grid_Yp_cm,Vtheta,15)%,'LineStyle','none')
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
Vtheta_fft=fft(Vtheta,[],1)/Nt; % --> Vtheta_fft(m,r)
%-----------Vr
Vr_fft=fft(Vr,[],1)/Nt; % --> Vtheta_fft(m,r)
%-----------modes
Nm = Nt;
M=ceil((Nm+1)/2);
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
%  ###########################################################################################################################
%  ################################################# End Time Loop ###########################################################
%  ###########################################################################################################################





%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   RECOMPOSITION TO SPATIAL DOMAIN AND PLOTS
%  ###########################################################################################################################
%  ###########################################################################################################################
% From here we first plot the energetic and theoretical spectra using dimensions 
% of the experiments. We also rebuid the velocity field to evaluate the
% errors between Besse-Fourier decomposition and the initial field. The
% integral of the energy (i.e 1/S * int_0^R int_0^(2*pi/alpha) 0.5*(u^2 + v^2) r dr dtheta, with S the surface)
% is also perfomed in the spectral and physical domain to be sure that we 
% retrieve the same quantities.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                SPECTRA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -------------------------------- Theoretical spectra
Cz = 0.5;
Ck = 6;
beta = beta{cexp,2}/100.; % cm-1.s-1
epsilon = 1*10.^(-5); % cm-2.s-3
EZ_Theo = (Cz*beta.^2.).*(1./R).*(J_root(1,1:Nrk)./R).^(-5.);
ER_Theo = Ck.*(epsilon.^(2./3.))*(1./R).*(J_root(2,1:Nrk)./R).^(-5./3.);
n = [1:Nrk]; % non-dimensional
U=0.0044*100;%cm
L_beta=(epsilon./beta.^3.).^(1./5.);
Lo_beta = (Ck/Cz).^(3./10.).*(J_root(1,1:Nrk)./J_root(2,1:Nrk)).^(1./2.);
Lhat_beta = (J_root(1,1:Nrk)./n).*Lo_beta.*L_beta;
L_R = sqrt(2.*U./beta); 
n_r=R/L_R;
% -------------------------------- plots spectra
figure;hold on
plot(n,mean(EZn,2),'r','linewidth',2.)
hold on
plot(n,mean(ERn,2),'k','linewidth',2.)
plot(n,EZ_Theo,'--r')
plot(n,ER_Theo,'--k')
% line([10 10], get(gca, 'ylim'));
set(gca,'YScale', 'log','XScale', 'log')
ylim([10.^(-7.) 0.1])
xlim([n(1) n(end)])
box on
xlabel('n','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E (cm^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
set(gca,'FontSize',13,'FontWeight','bold')
title(['\epsilon = ',num2str(epsilon),' cm^2.s^{-3}'],'FontSize',11)



EnergyT = sum(mean(EZn,2)) + sum(mean(ERn,2));
disp('####################################')
disp(' Energy from spectral coefficients')
disp('Energy = '), disp(num2str(EnergyT))

%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   PROJECTION IN SPACIAL DOMAIN
%  ###########################################################################################################################
%  ###########################################################################################################################

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             Normalisation Bessel error:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% questa mapa permette di verificare che la normalisazione dei funzioni di
% Bessel sia sotto un certo errore. Per minimizare quello errore Norm_mn
% deve essere il piu vicino possibile a 1. Questa mapa permette di
% scegliere il numero di Nrk che condusce a un erorre decente.
figure; hold on
contourf(Norm_mn)
title(['Normalisation errors = ',num2str(rms(rms(Norm_mn)))])
colorbar

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                   PROJECTION ON BESSEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vtheta_fft_proj=zeros(Nm,Nr);
Vr_fft_proj=zeros(Nm,Nr);
% for positive frequences m
for m=1:M%Nt2
for ir=1:Nr
Vtheta_fft_proj(m,ir) = sum( Cmn_U(m,:) .* besselj(m-1,J_root(m,1:Nrk).*r(ir)/R));
Vr_fft_proj(m,ir) = sum( Cmn_V(m,:) .* besselj(m-1,J_root(m,1:Nrk).*r(ir)/R));
end
end
% for negative frequences m
mm=0;
for m=Nm:-1:M+1
   mm=mm+1;
for ir=1:Nr
Vtheta_fft_proj(m,ir) = sum( Cmn_U(m,:) .* besselj(-mm,J_root(mm+1,1:Nrk).*r(ir)/R));
Vr_fft_proj(m,ir) = sum( Cmn_V(m,:) .* besselj(-mm,J_root(mm+1,1:Nrk).*r(ir)/R));
end 
end
% -------------------------------- plots
figure; hold on
plot(real(Vtheta_fft(1,:)))
plot(real(Vtheta_fft_proj(1,:)),'*')
title('Reconstruct: one m=0 along radius')
figure; hold on
plot(real(Vtheta_fft(10,:)))
plot(real(Vtheta_fft_proj(10,:)),'*')
title(['Reconstruct: one m=9 along radius, rms(error) = ',num2str(rms(rms(real(Vtheta_fft)-real(Vtheta_fft_proj))))])
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
alpha=round((2*pi)/(theta(end)-theta(1)));
disp(['Attenzione m = ',num2str(alpha),'*m'])
%
Vtheta_proj= zeros(size(Vtheta));
Vtheta_proj_BF= zeros(size(Vtheta));
Vr_proj_BF= zeros(size(Vr));
% loop sui raggi e theta per rimappare
for ir=1:Nr
    for itheta=1:Nt
       % Vtheta_proj(itheta,ir)= sum(Vtheta_fft(:,ir) .* exp((-2.*pi*i)/M.*-(m-1).*(itheta-1))); 
        Vtheta_proj(itheta,ir)= sum(Vtheta_fft(:,ir) .* exp((1i).*(alpha.*m).*theta(itheta))); 
        Vtheta_proj_BF(itheta,ir)= sum(Vtheta_fft_proj(:,ir) .* exp((1i).*(alpha.*m).*theta(itheta))); 
        Vr_proj_BF(itheta,ir)= sum(Vr_fft_proj(:,ir) .* exp((1i).*(alpha.*m).*theta(itheta))); 
    end
end
% -------------------------------- plots maps
if(0==0)
figure; hold on
contourf(Grid_Xp_cm,Grid_Yp_cm,Vtheta,15)%,'LineStyle','none'
title('Initial zonal velocity')
colorbar
% figure; hold on
% contourf(Grid_Xp_cm,Grid_Yp_cm,real(Vtheta_proj))
% title('Zonal velocity from Fourier coeffs')
figure; hold on
contourf(Grid_Xp_cm,Grid_Yp_cm,real(Vtheta_proj_BF),15)
colorbar
error_fit = rms(rms(real(Vtheta)-real(Vtheta_proj_BF)));
title(['Zonal velocity from Fourier-Bessel coeffs, rms(error) = ',num2str(error_fit),'m.s-1'])
end
% -------------------------------- plots profils
figure
plot(real(Vtheta(:,100)))
hold on
plot(real(Vtheta_proj(:,100)),'*')
plot(real(Vtheta_proj_BF(:,100)),'-*')
title('Zonal Velocity profils along azimuth')

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
disp(' Energy from physical velocity')
disp('Energy = '), disp(num2str(IAnalytique))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for itheta=1:Nt
    IA_rT = r.* 0.5.* (abs(Vtheta_proj_BF(itheta,:)).^2 + abs(Vr_proj_BF(itheta,:)).^2);
    IA_T(itheta) = trapz(r,IA_rT);
end
IA = trapz(theta,IA_T)./A;

% I = trapz(y,trapz(x,F,2))
   
disp('####################################')
disp(' Energy from rebuilt physical velocity')
disp('Energy = '), disp(num2str(IA))



%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   SAVE
%  ###########################################################################################################################
%  ###########################################################################################################################
if(ifsave==1)
fileout1 = [roots,names{cexp,2},'/',output,'/Ezn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(nFrames),'_window',num2str(windowing)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(EZn,1),'ulong');
fwrite(fid,size(EZn,2),'ulong');
fwrite(fid,EZn(:),'float');
fclose(fid);
% % % 
fileout2 = [roots,names{cexp,2},'/',output,'/ERn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(nFrames),'_window',num2str(windowing)];
filename2 = sprintf('%s.mtx',fileout2);
fid = fopen(filename2,'wb');
fwrite(fid,size(ERn,1),'ulong');
fwrite(fid,size(ERn,2),'ulong');
fwrite(fid,ERn(:),'float');
fclose(fid);
% % % 
R_trunc=R;
save([roots,names{cexp,2},'/',output,'/E_infos.mat'],'n','epsilon','Ck','Cz','L_beta','R_trunc','error_fit','Nrk','M','r','theta','Nrmin','Ntmin','Nrmax','Ntmax','Nr','Nt')
end