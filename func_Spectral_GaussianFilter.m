function [UfilterFG,VfilterFG,EG,kx,ky,Nkx,Nky] = func_Spectral_GaussianFilter(l,Nx,Ny,U,V,gaussian)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          FOURIER DECOMPOSITION for U:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In Matlab, the FFT is not normalized by the number of point Nx and Ny of
% the initial velocity field. The IFFT is normalizated. Then, if you
% perfrom only FFT it has to be normalizated if you do back and forth then
% you don't have to. The number of modes is Nkx=ceil((Nx+1)/2), with Nx and
% Ny being odd or even.
%-----------FT(U) sur x
U_fft=fft(U,[],2); % !! I do not normalize
V_fft=fft(V,[],2); % !! I do not normalize
Nkx=ceil((Nx+1)/2);
%-----------FT(U) sur y
U_2fft=fft(U_fft,[],1); % !! I do not normalize
V_2fft=fft(V_fft,[],1); % !! I do not normalize
Nky=ceil((Ny+1)/2);
%----------- Vector of modes ----------------------------------------------
% The vector of modes kx = [0...Nkx-1, Fq of Lunquist, -(Nkx-1)...-1] if Nx 
% and Ny are even and kx = [0...Nkx-1, -(Nkx-1)...-1] if odd.
%--------------------------------------------------------------------------
% vecteur des modes kx:
kx = zeros(1,Nx);
kx(1,1:Nkx) = [1:Nkx]-1;
if(mod(Nx,2)==0) % even
kx(1,Nkx+1:Nx)=-(Nkx-2:-1:1);
else % odd
kx(1,Nkx+1:Nx)=-(Nkx-1:-1:1);
end
% vecteur des modes ky:
ky = zeros(Ny,1);
ky(1:Nky) = [1:Nky]-1;
if(mod(Ny,2)==0) % even
ky(Nky+1:Ny)=-(Nky-2:-1:1);
else % odd
ky(Nky+1:Ny)=-(Nky-1:-1:1);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                      Filtrage spectral:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The spectral filtering is the product of the spectral field, obtains in 
% the preceding spectral decomposition of the velocity field, with the
% gaussian filter decomposed in the spectral space. To compute another
% function than the gaussian proposed below one has to find its fourier
% transform FT as in the exemple shown below and enter its expression for
% TFG = .... and TFG_toNorm that differ by the terms (Mkx.*Nx).^2) + ((Mky.*Ny).^2). 
% Here is the analytical expression of a guassian G and its fourier transform
% FT(G). Used as a filter the product of the spectral modes
% with the spectral gaussian leads to a gaussian filter in the same way than
% a convolution product of a gaussian with a velocity field. Spectral
% filtering is instead way much faster.
% The spectral gaussian is 2D and its expression is: 
% G     = (6/pi)^1/2 /l^2 * exp(-6 * (gx^2 + gy^2) / l.^2 );
% FT(G) = (6/pi)^1/2 /l^2 * 2pi/A*A * exp( -l^2 (kx^2 + ky^2)/24
% with A = 2*3^1/2*(1/l^2)^1/2. This expression is easy to find by
% computing the fft of G in wolfram alpha. Note that kx and ky are 2D
% matrices. Here l is the cutoff frequency and is analogue at the gaussian
% standard deviation sigma = l/(12)^1/2. Note that other gaussian function
% are avalaible.  
% !!! Within the first length scale l=1,2,3,.. the normalisation is
% underestimated and gives weird results.
Mkx = ones([Nx]).*kx.*2.*pi./Nx; % In the fourier transform the modes k=[0....pi]
Mky = ones([Ny]).*ky.*2.*pi./Ny; % In the fourier transform the modes k=[0....pi]
%--------------------------------------------------------------------------------
%--------------------- Gaussian of the type exp(-6r^2/l^2) from Physical mechanism 
% of the inverse energy cascade of two-dimensional turbulence: a numerical investigation 
% Z. XIAO et al (2009), doi:10.1017/S0022112008004266
if(gaussian==1)
A = 2.*sqrt(3.)*sqrt(1./(l.^2.));
C = sqrt(6./pi)./(l.^2);
TFG = C.*((2.*pi)./(A.*A)).*exp(-( (l.^2.).* ((Mkx.^2) + (Mky.^2) )./24));
TFG_toNorm = C.*((2.*pi)./(A.*A)).*exp(-( (l.^2.).* (((Mkx.*Nx).^2) + ((Mky.*Ny).^2) )./24));
%--------------------------------------------------------------------------------
%--------------------- Gaussian of the type exp(-9r^2/(2l^2)) Altimetric measurements 
% of rotating turbulence: cyclone-anticyclone asymmetry, inertial and Kelvin waves 
% and spectral characteristics, Y. D. Afanasyev (2013), [http://dx.doi.org/10.1063/1.4826477]
elseif(gaussian==2)
A = 3.*sqrt(1./(l.^2.));
C = 9./(2.*pi.*l.^2);
TFG = C.*((2.*pi)./(A.*A)).*exp(-( (l.^2.).* ((Mkx.^2) + (Mky.^2) )./18));
TFG_toNorm = C.*((2.*pi)./(A.*A)).*exp(-( (l.^2.).* (((Mkx.*Nx).^2) + ((Mky.*Ny).^2) )./18));
end
%--------------------------------------------------------------------------------
% In the exponential k=[0....pi] and cannot be in pixels while the integral
% of the gaussian to normalize is the integrals over all modes which then
% set the total energy in the gaussian to be unity. This integral is, i.e.
% sum{ sum{ TFG_Norm }} is constant has the integral of the gaussian has to
% be the same for all scale l, it means that the guassian goes from a dirac
% to well spread function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0==1)
figure
% imagesc(G./sum(sum(G)))
imagesc(TFG./sum(sum(TFG_toNorm)))
colorbar
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UG_2fft = U_2fft.*TFG./(sum(sum(TFG_toNorm))); % Ici sum{sum{g}} est la Norm du filtre gaussien G pour que la guassienne soit bornee a 1
VG_2fft = V_2fft.*TFG./(sum(sum(TFG_toNorm))); % Ici sum{sum{g}} est la Norm du filtre gaussien G pour que la guassienne soit bornee a 1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Energy:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Here we compute the energy of the filtred spectra to check how the cut
% off frequency act on the energetic distribution. The correponding modes
% ar kx and ky in pixels and l = 2pi/k, in pixels as well. 
% E_fil = 0.5.* (U_fil_2fft.*conj(U_fil_2fft) + V_fil_2fft.*conj(V_fil_2fft));
EG = 0.5.*(UG_2fft./(Nx.*Ny)).*conj(UG_2fft./(Nx.*Ny)); % Ici Nx*Ny est la Normalisation de la FFT qui est deja integre a la proj plus bas
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                  PROJECTION ON FOURIER:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% La projection avec ifft est plus precise que celle obtenue en
% reconstruisant a la main avec des boucles sur exp(2*pi*i*k*x).
% Voir partie commente plus bas.
% ifft sur y
UG_ifft = ifft(UG_2fft,[],1); %% ifft normalizes by Ny
VG_ifft = ifft(VG_2fft,[],1); %% ifft normalizes by Ny
% ifft sur x
UfilterFG = ifft(UG_ifft,[],2); %% ifft normalizes by Nx
VfilterFG = ifft(VG_ifft,[],2); %% ifft normalizes by Nx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It goes much faster to do the ifft than using the following projection to
% spatial space by doing u(y,x) = FT(u) * e^2iTTkx. It however the same 
% calculation.
% % % % % % U_filterFG_ky = zeros(size(V));
% % % % % % U_filterFG= zeros(size(U));
% % % % % % % projection sur ky pour revenir en spatial sur y
% % % % % % for iy=1:Ny
% % % % % %     for ikx=1:Nx
% % % % % %         % -> Si fft non normalisee par N
% % % % % %         U_filterFG_ky(iy,ikx)= (1./Ny).*sum(UG_2fft(:,ikx) .* exp((2.*pi*1i)/Ny.*(ky).*(iy))); 
% % % % % %         % -> Si fft est normalisee par N
% % % % % % %         U_filterFG_ky(iy,ikx)= sum(UG_2fft(:,ikx) .* exp((2.*pi*1i)/Ny.*(ky).*(iy))); 
% % % % % %     end
% % % % % % end
% % % % % % % projection sur kx pour revenir en spatial sur x
% % % % % % for iy=1:Ny
% % % % % %     for ix=1:Nx
% % % % % %         % -> Si fft non normalisee par N
% % % % % %         U_filterFG(iy,ix)= (1./Nx).*sum(U_filterFG_ky(iy,:) .* exp((2.*pi*1i)/Nx.*(kx).*(ix)));
% % % % % %         % -> Si fft est normalisee par N
% % % % % % %         U_filterFG(iy,ix)= sum(U_filterFG_ky(iy,:) .* exp((2.*pi*1i)/Nx.*(kx).*(ix)));
% % % % % %     end
% % % % % % end
end