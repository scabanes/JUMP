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

%% INPUTS: 
% # Set PDS parameters...
% # 'Vtheta_tot' and 'Vr_tot' the 2D velocity fields of azimuthal and
%    radial velocity of size (Nti*Nri,time).
% # Nri and Nti radial and azumuthal point numbers before truncation.
% # Ntmin and Ntmax are truncated azimuthal indices.
% # Nrmin and Nrmax are truncated radial indices.

%% OUTPUTS:
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                   LOAD:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Carichiamo input--------------------------
[Vtheta_tot,b]=loadmtx([roots,Name,NameVt]);
[Vr_tot,b]=loadmtx([roots,Name,NameVr]);
[returnOK] = Maps(Vtheta_tot,Vr_tot,Grid_Xp_cm,Grid_Yp_cm,GridR_2C,GridT_2C,Nti,Nri,Nrmin,Nrmax,Ntmin,Ntmax,theta,1)
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
EZn = zeros(Nrk,nTime);%Tmax);
ERn = zeros(Nrk,nTime);%Tmax);
Emn_t = zeros(Nrk*Nt,nTime);%Tmax);
%  ###########################################################################################################################
%  ##################################################### Time Loop ###########################################################
%  ###########################################################################################################################
iit=0;
for it=Itime:Tmax
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
Vtheta=Vtheta(Ntmin:Ntmax,Nrmin:Nrmax).*Tocm;
Vr=Vr(Ntmin:Ntmax,Nrmin:Nrmax).*Tocm;
URMS(it) = rms(rms(Vtheta + Vr));
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
[Emn,Norm_mn,M,Nm,Cmn_U,Cmn_V,Vtheta_fft,Vr_fft] = FourierBesselDecomp(Vtheta,Vr,Nt,Nr,Nrk,r,R);

iit = iit+1;
EZn(:,iit) = Emn(1,:);
ERn(:,iit) = sum(Emn(2:end,:));
Emn_t(:,iit) = Emn(:);
 end
%  ###########################################################################################################################
%  ################################################# End Time Loop ###########################################################
%  ###########################################################################################################################



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
load('besselzeros2_C.mat'); 

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
plot(real(Vtheta(:,10)))
hold on
plot(real(Vtheta_proj(:,10)),'*')
plot(real(Vtheta_proj_BF(:,10)),'-*')
title('Zonal Velocity profils along azimuth')

%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   ENERGY INTEGRALS
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
E_Ispectral = sum(mean(EZn,2)) + sum(mean(ERn,2));
disp('------------------------------------')
Urms = sqrt(E_Ispectral);
disp('Urms = '), disp([num2str(Urms),' cm/s'])
disp('####################################')
disp(' Energy from spectral coefficients')
disp('Energy = '), disp(num2str(E_Ispectral))

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
fileout1 = [roots,Name,'/Ezn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nTime),'_Itime_',num2str(Itime)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(EZn,1),'ulong');
fwrite(fid,size(EZn,2),'ulong');
fwrite(fid,EZn(:),'float');
fclose(fid);
% % % 
fileout2 = [roots,Name,'/ERn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nTime),'_Itime_',num2str(Itime)];
filename2 = sprintf('%s.mtx',fileout2);
fid = fopen(filename2,'wb');
fwrite(fid,size(ERn,1),'ulong');
fwrite(fid,size(ERn,2),'ulong');
fwrite(fid,ERn(:),'float');
fclose(fid);
% % % 
fileout2 = [roots,Name,'/Emn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nTime),'_Itime_',num2str(Itime)];
filename2 = sprintf('%s.mtx',fileout2);
fid = fopen(filename2,'wb');
fwrite(fid,size(Emn_t,1),'ulong');
fwrite(fid,size(Emn_t,2),'ulong');
fwrite(fid,Emn_t(:),'float');
fclose(fid);
% % % 
save([roots,Name,'/SpectralAnalysis_infos.mat'],'Urms','R','Nrk','M','nFrames','r','theta','Nrmin','Ntmin','Nrmax','Ntmax','Nr','Nt')
end