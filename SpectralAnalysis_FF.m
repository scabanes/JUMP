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
% Locally cartesian grid see Read. 2015
GridR_2C=GridR_2C(Ntmin:Ntmax,Nrmin:Nrmax);
GridT_2C=GridT_2C(Ntmin:Ntmax,Nrmin:Nrmax);
% -----------------------------------------------------Truncation
theta = theta(1:length(Ntmin:Ntmax)); % Attenzione theta deve partire da zero..
% !!!!!! Attention encore un doute
r = r(Nrmin:Nrmax); % Attenzione theta deve partire da zero..
% R=r(end);
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
% % % % WA = Windowing.WindowingTukey(Grid_Xp_cm,1); % Tukey function on the Azimuthal direction.
% % % % WR = Windowing.WindowingTukey(Grid_Xp_cm,2); % Tukey function on the Radial direction.
%  ###########################################################################################################################
%  ##################################################### Time Loop ###########################################################
%  ###########################################################################################################################
iit=0;
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
Vtheta=Vtheta(Ntmin:Ntmax,Nrmin:Nrmax);%.*Tocm;
Vr=Vr(Ntmin:Ntmax,Nrmin:Nrmax);%.*Tocm;
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
% if(windowing==1) 
% % ------------------------Azimuthal windowing
% Vtheta = Vtheta.*WA;
% Vr = Vr.*WA;
% % ------------------------Radial windowing
% Vtheta = Vtheta.*WR;
% Vr = Vr.*WR;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(0==0)
% figure; title('section we are working on')
% contourf(GridT_2C,GridR_2C,Vtheta,250,'LineStyle','none')
% end
%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   DECOMPOSITION IN SPECTRAL DOMAIN
%  ###########################################################################################################################
%  ###########################################################################################################################
% [r,Nx,Ny] = FourierFourierDecomp(Vtheta,Vr,r,theta,Ro);
FourierFourierDecomp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SAVE in Time
iit = iit+1;
%  Zonostrophic formulation
EZ_ky_t(:,iit) = EZ_ky;
ER_ky_t(:,iit) = ER_ky;
% QNSE
E1_kx_t(:,iit) = E1_kx;
E1_ky_t(:,iit) = E1_ky;

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                  PROJECTION ON FOURIER:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_proj= zeros(size(Vtheta));
U_proj_ky = zeros(size(Vtheta));
% projection sur ky pour revenir en spatial sur y
for iy=1:Ny
    for ikx=1:Nx
        % -> Si fft non normalisee par N
%         U_proj_ky(iy,ikx)= (1./Ny).*sum(U_2fft(:,ikx) .* exp((2.*pi*1i)/Ny.*(ky).*(iy))); 
        % -> Si fft est normalisee par N
        U_proj_ky(ikx,iy)= sum(U_2fft(ikx,:) .* exp((2.*pi*1i)/Ny.*(ky).*(iy))); 
    end
end
% projection sur kx pour revenir en spatial sur x
for iy=1:Ny
    for ix=1:Nx
        % -> Si fft non normalisee par N
%         U_proj(iy,ix)= (1./Nx).*sum(U_proj_ky(iy,:) .* exp((2.*pi*1i)/Nx.*(kx).*(ix)));
        % -> Si fft est normalisee par N
        U_proj(ix,iy)= sum(U_proj_ky(:,iy) .* exp((2.*pi*1i)/Nx.*(kx).*(ix)));
    end
end
if(0==1)
 figure; title('projected map')
 contourf(GridT_2C,GridR_2C,flip(real(U_proj)), 250, 'LineStyle','none');
 colorbar 
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                      SURFACE INTEGRALS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------- from spectral data
% -> si la fft n est pas normalise alors pour l energie somme dans les deux
% direction dx et dy on a un facteur 1/((Nx*Ny).^2)
% IS = 0.5.*(1/((Nx*Ny).^2)).* sum(sum(U_2fft.*conj(U_2fft)))
% -> lorsque la fft est normalisee par Nx et Ny on a pas ce facteur
IS = 0.5.* sum(sum(U_2fft.*conj(U_2fft)))

%------------- from analytical data
% premiere integrale optionel qui somme sur les points x et y et non dans l
% espace et donc divise par le nombre de point
IA = 0.5.*(1/((Nx*Ny))).*sum(sum(Vtheta.*Vtheta))
% Integrale sur l espace x et y et divise par la surface
S = (x(end)-x(1)).*(y(end)-y(1));
IAA = (1/S).*0.5.*trapz(y,trapz(x,Vtheta.*Vtheta,1))

%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   SAVE
%  ###########################################################################################################################
%  ###########################################################################################################################
if(ifsave==1)
fileout1 = [roots,Name,'/EZ_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(EZ_ky_t,1),'ulong');
fwrite(fid,size(EZ_ky_t,2),'ulong');
fwrite(fid,EZ_ky_t(:),'float');
fclose(fid);
% % % 
fileout1 = [roots,Name,'/ER_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(ER_ky_t,1),'ulong');
fwrite(fid,size(ER_ky_t,2),'ulong');
fwrite(fid,ER_ky_t(:),'float');
fclose(fid);
% % % 
fileout1 = [roots,Name,'/E1_kx_Fr_',num2str(nTime),'_window_',num2str(windowing)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(E1_kx_t,1),'ulong');
fwrite(fid,size(E1_kx_t,2),'ulong');
fwrite(fid,E1_kx_t(:),'float');
fclose(fid);
% % % 
fileout1 = [roots,Name,'/E1_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(E1_ky_t,1),'ulong');
fwrite(fid,size(E1_ky_t,2),'ulong');
fwrite(fid,E1_ky_t(:),'float');
fclose(fid);
% % % 
save([roots,Name,'/SpectralAnalysis_FF_infos.mat'],'Nkx','Nky','nFrames','x','y')
end