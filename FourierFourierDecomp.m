% function [r,Nx,Ny] = FourierFourierDecomp(Vtheta,Vr,r,theta,Ro)
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
x = Ro.*theta; 
y = r;
Nx = length(x);
Ny = length(y);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                    FOURIER DECOMPOSITION ALONG X and Y:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% la fft deve essere normalisata per il numero di punti in azimuth che corrisponde
% al numero dei modi Nm che escono della fft. Tra questi Nm modi ci sono M modi
% zonali m che corrisponde al numero dei punti unici simetrici per M ne l vettore _fft.
%-----------U
U_fft=fft(Vtheta,[],1)/Nx; % !! it is normalized
V_fft=fft(Vr,[],1)/Nx; % !! it is normalized
Nkx=ceil((Nx+1)/2);
%
U_2fft=fft(U_fft,[],2)/Ny; % !! it is normalized
V_2fft=fft(V_fft,[],2)/Ny; % !! it is normalized
Nky=ceil((Ny+1)/2);
% vecteur des modes kx
kx = zeros(Nx,1);
kx(1:Nkx) = [1:Nkx]-1;
if(mod(Nx,2)==0) % even
kx(Nkx+1:Nx)=-(Nkx-2:-1:1);
else % odd
kx(Nkx+1:Nx)=-(Nkx-1:-1:1);
end
% vecteur des modes ky
ky = zeros(1,Ny);
ky(1,1:Nky) = [1:Nky]-1;
if(mod(Ny,2)==0) % even
ky(1,Nky+1:Ny)=-(Nky-2:-1:1);
else % odd
ky(1,Nky+1:Ny)=-(Nky-1:-1:1);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Energy:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pour l energie on a 0.5 * U^2
E = 0.5.* (U_2fft.*conj(U_2fft) + V_2fft.*conj(V_2fft)); %E(x,y)
E1 = 0.5.* (U_2fft.*conj(U_2fft)); %E(x,y)
E2 = 0.5.* (V_2fft.*conj(V_2fft)); %E(x,y)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Energy pour ky et -ky:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EZ_Tky = E(1,:);
ER_Tky = sum(E(2:end,:),1);
%
E1_Tkx = sum(E1(:,:),2);
E1_Tky = sum(E1(:,:),1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Energy pour les ky:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Je somme ici les modes ky et -ky qui contribue tout deux a l ernergie du 
% mode ky
EZ_ky = zeros([Nky,1]);
ER_ky = zeros([Nky,1]);
%
E1_ky = zeros([Nky,1]);
% ce vecteur est juste pour organiser la collecte de modes.
vky = ky.*sign(ky); 
for iky_p=0:Nky-1 % que les ky positifs
    l=(vky==iky_p);
    EZ_ky(iky_p+1) = sum(EZ_Tky(l)); % j'aurai pu ajouter simplement un *2
    ER_ky(iky_p+1) = sum(ER_Tky(l)); % j'aurai pu ajouter simplement un *2
    %
    E1_ky(iky_p+1) = sum(E1_Tky(l)); % j'aurai pu ajouter simplement un *2
end
%
E1_kx = zeros([Nkx,1]);
% E2_kx = zeros([Nkx,1]);
% ce vecteur est juste pour organiser la collecte de modes.
vkx = kx.*sign(kx); 
for ikx_p=0:Nkx-1 % que les ky positifs
    l=(vkx==ikx_p);
    E1_kx(ikx_p+1) = sum(E1_Tkx(l)); % j'aurai pu ajouter simplement un *2
%     E2_kx(ikx_p+1) = sum(E2_Tkx(l)); % j'aurai pu ajouter simplement un *2
end

% end

