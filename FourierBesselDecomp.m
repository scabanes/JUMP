function [Emn,Norm_mn,M,Nm,Cmn_U,Cmn_V,Vtheta_fft,Vr_fft] = FourierBesselDecomp(Vtheta,Vr,Nt,Nr,Nrk,r,R)
%% INPUTS:
%-# Vtheta & Vr are azimuthal and radial 2D velocity fields.
%   Matrice [Vtheta] & [Vr] have the size Nt * Nr
%-# Nt is the number of points in azimuth theta
%-# Nr is the number of points in radius r
%-# Nrk is the number of radial wavenumber wanted in radius
%-# r is a vector of length Nr corresponding tho the radius (possibly
%   truncated)
%-# R is the maximum radius after radial truncation (if truncation is
%   applied)
%% OUTPUTS:
%-# Emn is a matrix containing the energy of spectral modes m & n
%   The size of [Emn] is Nm * Nrk
%-# Norm_mn is a normalisation matrix for spectral modes m & n
%   The size of [Norm_mn] is Nm * Nrk. This matrix should be unity if
%   normalisation is good.
%-# Cmn_U & Cmn_V are decomposition coefficients of spectral modes m & n
%   The size of [Cmn_X] is Nm * Nrk
%-# M is a vector of the azimuthal modes m
%-# Nm are the number of azimuthal modes
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
% ---- First we load the zeros of the bessel functions.
load('besselzeros2_C.mat'); 
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

end

