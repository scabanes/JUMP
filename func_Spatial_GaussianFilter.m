function [U_fil,V_fil,E_fil,kx,ky,Nkx,Nky] = func_Spatial_GaussianFilter(l,Nx,Ny,U,V,gaussian)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                           Filtering In space:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The spatial filter results from the convolution product of the 2D
% velocity field with a gaussian. The gaussian has to be normalized by the
% integral over both x and y direction in order to get unity.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          FILTERING In SPACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r1=Nx/2; % This size have have to be enlarge otherwise the gaussian is poorly resolved when filtering at large scale!
r2=Ny/2; % This size have have to be enlarge otherwise the gaussian is poorly resolved when filtering at large scale!
gx = [-r1:r2];
gy = [-r1:r2]';
Ngx=length(gx);
Ngy=length(gy);
% Gaussienne I.
if(gaussian==1)
C = sqrt(6./pi)./(l.^2);
g = C.*exp(-6.*(gx.^2 + gy.^2)/(l.^2));
% Gaussienne II.
elseif(gaussian==2)
C = 9./(2.*pi.*(l.^2));
g = C.*exp(-9.*(gx.^2 + gy.^2)/(2.*l.^2));
end
% Ici sum{sum{g}} est la Norm du filtre gaussien G pour que la guassienne 
% soit bornee a 1
hg = g./sum(sum(g));
% produit de convolution: filtrage gaussien
U_fil = conv2(U,hg,'same');
V_fil = conv2(V,hg,'same');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                          Spectral decomposition of a
%                                          guassian with FFT to check the
%                                          analysical expression used in 
%                                          func_Spectral_GaussianFilter.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With the discrete fourier transform of the gaussian, we retrieve the
% estimat TFG used in func_Spectral_GaussianFilter.m, the abs(hg_2fft) maps
% can be directly compared with TFG. Only for l<3 or 4 this is not true!!
% Fourier decomposition de la gaussienne
hg_fft = fft(hg,[],2);%/(r1*2+1);
Nkgx=ceil((Ngx+1)/2);
hg_2fft = fft(hg_fft,[],1);%/(r2*2+1);
Nkgy=ceil((Ngy+1)/2);
% vecteur des modes kx
kgx = zeros(1,Ngx);
kgx(1,1:Nkgx) = [1:Nkgx]-1;
if(mod(Ngx,2)==0) % pair
kgx(1,Nkgx+1:Ngx)=-(Nkgx-2:-1:1);
else % impair
kgx(1,Nkgx+1:Ngx)=-(Nkgx-1:-1:1);
end
% vecteur des modes ky
kgy = zeros(Ngy,1);
kgy(1:Nkgy) = [1:Nkgy]-1;
if(mod(Ngy,2)==0) % pair
kgy(Nkgy+1:Ngy)=-(Nkgy-2:-1:1);
else % impair
kgy(Nkgy+1:Ngy)=-(Nkgy-1:-1:1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0==1)
figure
imagesc(abs(hg_2fft))
colorbar
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       FOURIER DECOMPOSITION  for U_fil:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we estimate E_fil to obtain the spectrum of the velocity field once
% it has been filterd.
% la fft deve essere normalisata per il numero di punti in azimuth che corrisponde
% al numero dei modi Nm che escono della fft. Tra questi Nm modi ci sono M modi
% zonali m che corrisponde al numero dei punti unici simetrici per M ne l vettore _fft.
%-----------U
U_fil_fft=fft(U_fil,[],2)/Nx; % !! it is normalized
V_fil_fft=fft(V_fil,[],2)/Nx; % !! it is normalized
Nkx=ceil((Nx+1)/2);
%
U_fil_2fft=fft(U_fil_fft,[],1)/Ny; % !! it is normalized
V_fil_2fft=fft(V_fil_fft,[],1)/Ny; % !! it is normalized
Nky=ceil((Ny+1)/2);
%----------- Vector of modes ----------------------------------------------
% The vector of modes kx = [0...Nkx-1, Fq of Lunquist, -(Nkx-1)...-1] if Nx 
% and Ny are even and kx = [0...Nkx-1, -(Nkx-1)...-1] if odd.
%--------------------------------------------------------------------------
kx = zeros(1,Nx);
kx(1,1:Nkx) = [1:Nkx]-1;
if(mod(Nx,2)==0) % pair
kx(1,Nkx+1:Nx)=-(Nkx-2:-1:1);
else % impair
kx(1,Nkx+1:Nx)=-(Nkx-1:-1:1);
end
% vecteur des modes ky
ky = zeros(Ny,1);
ky(1:Nky) = [1:Nky]-1;
if(mod(Nx,2)==0) % pair
ky(Nky+1:Ny)=-(Nky-2:-1:1);
else % impair
ky(Nky+1:Ny)=-(Nky-1:-1:1);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Energy:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pour l energie on a 0.5 * U^2
% E_fil = 0.5.* (U_fil_2fft.*conj(U_fil_2fft) + V_fil_2fft.*conj(V_fil_2fft));
E_fil = 0.5.*(U_fil_2fft.*conj(U_fil_2fft));

end