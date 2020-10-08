%% INPUTS: 
% # Set PDS parameters...
% # 'Vtheta_tot' and 'Vr_tot' the 2D velocity fields of azimuthal and
%    radial velocity of size (Nti*Nri,time).
% # Nri and Nti radial and azimuthal point numbers before truncation.
% # Ntmin and Ntmax are truncated azimuthal indices.
% # Nrmin and Nrmax are truncated radial indices.
% # x and y are cartesian coordinates
% # r and theta are polar coordinates

%% OUTPUTS:
% # E_ky_t and E_ky_t are kinetic energy along the wavenumbers ky and kx as 
%   a function of time t. They are matrices of E_ky_t(Nky,t) and E_kx_t(Nkx,t)
% # Nkx and Nky, are the number of modes in x and y
% # kx and ky, are non-dimantional modes that leads to wavenumbers in m-1
%   following 2*pi*kx/Lx, with Lx is the domain's length in x and in meter.
%   Same for ky in y.

% _________________________________________________________________________
it=0;
%----------------------- Torino
FieldFrom_Torino
%----------------------- 2DJupiter_Obs
% FieldFrom_2DJupiter_Obs
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
%
disp(['step = ',num2str(it),' on ',num2str(Tmax)])

%##########################################################################
%                                                                   Torino:
%##########################################################################
file = [roots,Name,'/merged.civ2/fig_',num2str(numt),'.nc'];
if isfile(file)%------------------------------------------------------------------------------------- Init: file exists
% compter les images
numt = numt+gapt;% File exists.
%----------------------- Torino
FieldFrom_Torino
%##########################################################################
%                                                             2DJupiter_Obs:
%##########################################################################
% file = [roots,Name,'/StatisticalData.nc'];
% if isfile(file)%------------------------------------------------------------------------------------- Init: file exists
% %----------------------- 2DJupiter_Obs
% FieldFrom_2DJupiter_Obs
%##########################################################################
%##########################################################################
%##########################################################################


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------- Truncation of the velocity fields
Vtheta = Vtheta(idymin:idymax,idxmin:idxmax)';
Vr = Vr(idymin:idymax,idxmin:idxmax)';
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %% % % % end%------------------------------------------------ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if(0==0)
% figure; title('section we are working on')
% contourf(Grid_Xp_cm,Grid_Yp_cm,Vtheta,15)%,'LineStyle','none')
% colorbar
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if some data are missing in the velocity field we put a zero value
% instead of a NaN.
Vtheta(isnan(Vtheta))=0.;
Vr(isnan(Vr))=0.;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              WINDOWING:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The windowing make the fields to be periodic. Not necessarly used.
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
% Total energy
E_ky_t(:,iit) = E_ky;
%  Zonostrophic formulation
EZ_ky_t(:,iit) = EZ_ky;
ER_ky_t(:,iit) = ER_ky;
% QNSE
Ex_kx_t(:,iit) = Ex_kx;
Ey_kx_t(:,iit) = Ey_kx;
%
Ex_ky_t(:,iit) = Ex_ky;
Ey_ky_t(:,iit) = Ey_ky;

else% File does not exist.
    numt = numt + gapt;
end%------------------------------------------------------------------------------------- End: file exists
end%--------------------------------------------------------------------------------------End: LOOP on t

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
ETotal_S = 0.5.* sum(sum(U_2fft.*conj(U_2fft))) + 0.5.* sum(sum(V_2fft.*conj(V_2fft)))
%------------- from analytical data
% premiere integrale optionel qui somme sur les points x et y et non dans l
% espace et donc divise par le nombre de point
ETotal_A = 0.5.*(1/((Nx*Ny))).*sum(sum(Vtheta.*Vtheta)) + 0.5.*(1/((Nx*Ny))).*sum(sum(Vr.*Vr))
% Integrale sur l espace x et y et divise par la surface
S = (x(end)-x(1)).*(y(end)-y(1));
ETotal_AA = (1/S).*0.5.*trapz(y,trapz(x,Vtheta.*Vtheta,1)) + (1/S).*0.5.*trapz(y,trapz(x,Vr.*Vr,1))

disp('Energy integrate')
% ET = 0.5.*sum(sum((U_2fft.*conj(U_2fft) + V_2fft.*conj(V_2fft) ) ))
% ET = sum(sum(0.5.*U_2fft.*conj(U_2fft))) + sum(sum(0.5.*V_2fft.*conj(V_2fft)))
% ET = sum(EZ_ky_t(:,end),1) + sum(ER_ky_t(:,end),1)
disp(['EZ = ',num2str( mean(sum(EZ_ky_t,1)) )])
disp(['ER = ',num2str( mean(sum(ER_ky_t,1)) )])
disp(['ET = ',num2str( mean(sum(EZ_ky_t,1)) + mean(sum(ER_ky_t,1)) )])
disp(['ET_xy = ',num2str( mean(sum(Ex_ky,1)) + mean(sum(Ey_ky,1)) )])
%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   SAVE
%  ###########################################################################################################################
%  ###########################################################################################################################
if(ifsave==1)
fileout1 = [roots,Name,'/E_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(E_ky_t,1),'ulong');
fwrite(fid,size(E_ky_t,2),'ulong');
fwrite(fid,E_ky_t(:),'float');
fclose(fid);
%
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
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileout1 = [roots,Name,'/Ex_kx_Fr_',num2str(nTime),'_window_',num2str(windowing)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(Ex_kx_t,1),'ulong');
fwrite(fid,size(Ex_kx_t,2),'ulong');
fwrite(fid,Ex_kx_t(:),'float');
fclose(fid);
% % % 
fileout1 = [roots,Name,'/Ey_kx_Fr_',num2str(nTime),'_window_',num2str(windowing)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(Ey_kx_t,1),'ulong');
fwrite(fid,size(Ey_kx_t,2),'ulong');
fwrite(fid,Ey_kx_t(:),'float');
fclose(fid);
% % % %%%%%%%%%%%%%%%%%%%
fileout1 = [roots,Name,'/Ex_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(Ex_ky_t,1),'ulong');
fwrite(fid,size(Ex_ky_t,2),'ulong');
fwrite(fid,Ex_ky_t(:),'float');
fclose(fid);
% % % 
fileout1 = [roots,Name,'/Ey_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(Ey_ky_t,1),'ulong');
fwrite(fid,size(Ey_ky_t,2),'ulong');
fwrite(fid,Ey_ky_t(:),'float');
fclose(fid);
% % % 
save([roots,Name,'/SpectralAnalysis_FF_infos.mat'],'Nkx','Nky','nFrames','x','y','kx','ky','Ro')
end