%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   DECOMPOSITION IN SPECTRAL DOMAIN
%  ###########################################################################################################################
%  ###########################################################################################################################
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
E = 0.5 .* (U_2fft.*conj(U_2fft) + V_2fft.*conj(V_2fft)); %E(x,y)
E1 = 0.5 .* (U_2fft.*conj(U_2fft)); %E(x,y)
E2 = 0.5 .* (V_2fft.*conj(V_2fft)); %E(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Energy for ky et -ky:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_Tky = sum(E(:,:),1);
%
EZ_Tky = E(1,:);
ER_Tky = sum(E(2:end,:),1);
%
E1_Tky = sum(E1(:,:),1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Energy for kx et -kx:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1_Tkx = sum(E1(:,:),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Energy for ky:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we sum on modes ky and -ky that are both contributing to total energy 
% on ky modes
E_ky = zeros([Nky,1]);
%
EZ_ky = zeros([Nky,1]);
ER_ky = zeros([Nky,1]);
%
E1_ky = zeros([Nky,1]);
% A vector that corresponds to the different modes in matlab FFT.
vky = ky.*sign(ky); 
for iky_p=0:Nky-1 % que les ky positifs
    l=(vky==iky_p);
    E_ky(iky_p+1) = sum(E_Tky(l)); % This is equivalent than having a *2
    %
    EZ_ky(iky_p+1) = sum(EZ_Tky(l)); % This is equivalent than having a *2
    ER_ky(iky_p+1) = sum(ER_Tky(l)); % This is equivalent than having a *2
    %
    E1_ky(iky_p+1) = sum(E1_Tky(l)); % This is equivalent than having a *2
end
%
E1_kx = zeros([Nkx,1]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Energy for kx:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A vector that corresponds to the different modes in matlab FFT.
vkx = kx.*sign(kx); 
for ikx_p=0:Nkx-1 % que les ky positifs
    l=(vkx==ikx_p);
    E1_kx(ikx_p+1) = sum(E1_Tkx(l)); % This is equivalent than having a *2
%     E2_kx(ikx_p+1) = sum(E2_Tkx(l)); % This is equivalent than having a *2
end