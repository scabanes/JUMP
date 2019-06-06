clear all
close all
%**********************************************************
%           PARAMETRI DA SCEGLIERE
%**********************************************************
roots = '/media/simon/simon/ESP_29/'; % Root path..
run([roots,'InfosFile.m'])
%**********************************************************
% Dimensioni
Create_Grid_pol_Galperin2015
dr=R/Nri;
r=(1:Nri)*dr;
theta = [0:(2.*pi)/Nti:2.*pi-(2.*pi)/Nti];
%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   SPECTRAL DECOMPOSITION
%  ###########################################################################################################################
%  ###########################################################################################################################

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                   LOAD:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------Campo di velocita
[Vtheta_tot,a] = loadmtx([roots,NameVt]);
[Vr_tot,a] = loadmtx([roots,NameVr]);
%
if(isnan(Tmax)==1)
    Tmax=a(2);
end
nFrames = Tmax-Itime+1;
% -----------------------------------------------------Truncation
theta = theta(1:length(Ntmin:Ntmax)); % Attenzione theta deve partire da zero..
% !!!!!! Attention encore un doute
r = r(1:length(Nrmin:Nrmax)); % Attenzione theta deve partire da zero..
R=r(end);
% r = r(Nrmin:Nrmax); % Attenzione theta deve partire da zero..
% !!!!!! -------------------------
Nt = length(theta)
Nr = length(r)
%
[Vr] = Maps_Galperin2015(Vtheta_tot,Vr_tot,Grid_Xp_cm,Grid_Yp_cm,Nti,Nri,Nrmin,Nrmax,Ntmin,Ntmax,theta,1);
%
EZn = zeros(Nrk,nFrames);
ERn = zeros(Nrk,nFrames);
iit=0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              TIME LOOP:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for it=Itime:Tmax
display([num2str(it),'/',num2str(Tmax)])

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                              FOURIER BESSEL DECOMPOSITION OF 2D FIELDS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Emn,Norm_mn,M,Nm,Cmn_U,Cmn_V,Vtheta_fft,Vr_fft] = FourierBesselDecomp(Vtheta,Vr,Nt,Nr,Nrk,r,R);

iit = iit+1;
EZn(:,iit) = Emn(1,:);
ERn(:,iit) = sum(Emn(2:end,:));
end





%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   PROJECTION IN SPACIAL DOMAIN
%  ###########################################################################################################################
%  ###########################################################################################################################

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                   PROJECTION ON BESSEL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('besselzeros2_C.mat'); 

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


%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   ENERGY INTEGRALS
%  ###########################################################################################################################
%  ###########################################################################################################################

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                SPECTRA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------ Integral spectral de l energie
E_Ispectral = sum(mean(EZn,2)) + sum(mean(ERn,2));
disp('########### Spectral #########################')
disp('Energy = '), disp(num2str(E_Ispectral))
Urms = sqrt(E_Ispectral);
disp('Urms = '), disp([num2str(Urms),' cm/s'])

%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   ENERGY INTEGRALS
%  ###########################################################################################################################
%  ###########################################################################################################################

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

    
disp('############ Analytique ########################')
disp('Energy = '), disp(num2str(IAnalytique))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for itheta=1:Nt
    IA_rT = r.* 0.5.* abs(Vtheta_proj_step2(itheta,:)).^2;
    IA_T(itheta) = trapz(r,IA_rT);
end
IA = trapz(theta,IA_T)./A;

% I = trapz(y,trapz(x,F,2))
   
disp('########### Reconstructed #########################')
disp('Energy = '), disp(num2str(IA))


%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   SAVE
%  ###########################################################################################################################
%  ###########################################################################################################################
if(ifsave==1)
fileout1 = [roots,'/Ezn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nFrames),'_window_',num2str(windowing)];
filename1 = sprintf('%s.mtx',fileout1);
fid = fopen(filename1,'wb');
fwrite(fid,size(EZn,1),'ulong');
fwrite(fid,size(EZn,2),'ulong');
fwrite(fid,EZn(:),'float');
fclose(fid);
% % % 
fileout2 = [roots,'/ERn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nFrames),'_window_',num2str(windowing)];
filename2 = sprintf('%s.mtx',fileout2);
fid = fopen(filename2,'wb');
fwrite(fid,size(ERn,1),'ulong');
fwrite(fid,size(ERn,2),'ulong');
fwrite(fid,ERn(:),'float');
fclose(fid);
% % % 
save([roots,'/SpectralAnalysis_infos.mat'],'Urms','R','Nrk','M','nFrames','r','theta','Nrmin','Ntmin','Nrmax','Ntmax','Nr','Nt')
end