% Script per la determinazione della derivata prima e derivata seconda della Vz in dr
%
%% Input: - Vz.mtx
%         - Prof_Vz_sector_th1_th2.mtx

%% Output: - Prof_d1_Vz_dr.mtx, mappe dei profili della derivata prima di Vz in
%%           funzione di r, matrice [Nraggi*Ncerchi, totaltime] 

         % - Prof_d2_Vz_dr.mtx, mappe dei profili della derivata seconda di Vz in
%%           funzione di r, matrice [Nraggi*Ncerchi, totaltime] 

         % - Figura in cui si plottano il profilo medio di Vz (med temp e spaziale nel
%            settore circolare considerato), la derivata prima e la derivata seconda
%% Parametri da scegliere: -
clear all
%**********************************************************
%           PARAMETRI DA SCEGLIERE
%**********************************************************
cexp = 2;
expt = {'5';'6';'7';'8';'9';'3'}
roots = '../Dati-Torino/';
output = 'outputs'
%**********************************************************
% ##################################################################################################################################################
%                                                                                                                                 DATI INFORMAZIONI:
% ##################################################################################################################################################
EXPT_infos
load([roots,names{cexp,2},'/',output,'/infos.mat'])
nFrames = ((endi{cexp,2}-starti{cexp,2})/step + 1);
%nFrames=2150;   %%%solo per exp 8 altrimenti istruzione precedente!!!!!!
% ##################################################################################################################################################
%                                                                                                                                       BETA EFFECT:
% ##################################################################################################################################################
g = 9.81; %accelerazione di gravità m/s
Lx = 5; %dimensione del dominio in m
Ly = 5;
H0 = H0_exp{cexp,2}
Omega = omega_exp{cexp,2}
% omega = 0.206% in rad.s-1   QUESTO E' QUELLO FATTO DA SIMON
% H = [1.:120.]/215.%0.56 % in metri
%per il paraboloide e per la topografia di fondo
%Omega = Omega /(2*pi)%conversione in s^-1 dobbiamo farla??
f0=2*Omega;
s = Omega^2/(2*g);
Hc = H0-(s*(Lx^2+Ly^2)/12);
h=Hc+(r.^2).*(s);
h_fondo=(MaxRadius-r)*tan(11.*(2.*pi)/360.);
H=h-h_fondo;  
% ##################################################################################################################################################
%                                                                                                                                              LOAD:
% ##################################################################################################################################################
% name_mtx = ['Vz09_540_180_2160' ];
name_mtx = [roots,names{cexp,2},'/',output,'/Vz_',num2str(Nraggi),'_',num2str(Ncerchi),'_',num2str(nFrames)];
[Vz_int,b] = loadmtx(name_mtx);
% name1_mtx = ['Vr09_540_160_2160' ];
name1_mtx = [roots,names{cexp,2},'/',output,'/Vr_',num2str(Nraggi),'_',num2str(Ncerchi),'_',num2str(nFrames)];
[Vr_int,b] = loadmtx(name1_mtx);

totaltime = b(2);
dth = 360/Nraggi;

%r = (0:dr:r_max_cm);

dVz_r_time = zeros (Nraggi*Ncerchi,totaltime);
dVr_r_time = zeros (Nraggi*Ncerchi,totaltime);
VortRelat_time = zeros (Nraggi*Ncerchi,totaltime);
VortPot_time = zeros (Nraggi*Ncerchi,totaltime);
VortPot_zm_time = zeros (Ncerchi,totaltime);

% dVz_2_r_time = zeros (Nraggi*Ncerchi,totaltime);

%%% calcolo della derivata prima e seconda---------------------------------

for t=1:totaltime
    clear Vz
    clear Vr
    disp (t);
    Vz = reshape (Vz_int(:,t),Nraggi,Ncerchi);
    Vr = reshape (Vr_int(:,t),Nraggi,Ncerchi);
        
  dVz = zeros(Nraggi,Ncerchi);
  dVr = zeros(Nraggi,Ncerchi);

 % dVz_2 = zeros(Nraggi,Ncerchi);
  
  % calcolo le differenze sulla V zonale
    for i=2:Ncerchi-1
    
       %%% dal punto r(2) uso le differenze finite centrate fino al
       %%% punto r(Ncerchi-1); per il punto r(1) dVz =0;
       % derivata prima
       dVz(:,i) = Rpoints(i+1)*Vz(:,i+1) - Rpoints(i-1)*Vz(:,i-1);
       
       % derivata seconda
     %   dVz_2(:,i) = Vz(:,i+1) - 2*Vz(:,i) + Vz(:,i-1);
    
    end
 
     % calcolo le differenze sulla V radiale
    for ii=2:Nraggi-1
    
       %%% dal punto r(2) uso le differenze finite centrate fino al
       %%% punto r(Ncerchi-1); per il punto r(1) dVz =0;
       
       % derivata prima
       dVr(ii,:) = Vr(ii+1,:) - Vr(ii-1,:);
       
       % derivata seconda
     %   dVz_2(:,i) = Vz(:,i+1) - 2*Vz(:,i) + Vz(:,i-1);
    
     end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%% VETTORE-matrice DVz/Dr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% derivata prima
  dVz_dr = dVz./(2*dr);
% derivata prima
  dVr_dphi = dVr./(2*dphi);
  
  for j=1:Nraggi
  aa = ~isnan(Vz(j,:));
  aa_min = min(find(aa==1));
  aa_max = max(find(aa==1));
  %find(aa==1)
    %%per il punto r(Ncerchi) uso le differenze finite all'indietro
  dVz_dr (j,aa_max) = (Rpoints(aa_max).*Vz(j,aa_max) - Rpoints(aa_max-1).*Vz(j,aa_max-1))./(dr);
    %%per il punto r(0) uso le differenze finite in avanti
  dVz_dr (j,aa_min) = (Rpoints(aa_min+1).*Vz(j,aa_min+1) - Rpoints(aa_min).*Vz(j,aa_min))./(dr);
  end

    for jj=1:Ncerchi
  aaa = ~isnan(Vr(:,jj));
  aaa_min = min(find(aaa==1));
  aaa_max = max(find(aaa==1));
  %find(aa==1)
    %%per il punto r(Ncerchi) uso le differenze finite all'indietro
  dVr_dphi (aaa_max,jj) = (Vr(aaa_max,jj) - Vr(aaa_max-1,jj))./(dphi);
    %%per il punto r(0) uso le differenze finite in avanti
  dVr_dphi (aaa_min,jj) = (Vr(aaa_min+1,jj) - Vr(aaa_min,jj))./(dphi);
    end

    VortRelat = zeros(Nraggi,Ncerchi);
    VortPot = zeros(Nraggi,Ncerchi);
    VortPot1 = zeros(Nraggi,Ncerchi);
    VortPot_zm = zeros(Ncerchi,1);
    for k=1:Ncerchi
         VortRelat(:,k) = (1./Rpoints(k))*(dVz_dr(:,k) - dVr_dphi(:,k));
         VortPot(:,k) = (1./H(k)).* ( VortRelat(:,k) + f0); %2.*omega );  (1./H(k))
         %VortPot1(:,k) = (1./H(k)).* ( VortRelat(dd,k) + f0); %2.*omega );  (1./H(k))
         cc = ~isnan(VortPot(:,k));
%           if (cc)
%               VortPot1(:,k) = VortPot(cc,k)  %(1./H(k))
%           else 
%               VortPot1(:,k) = 0;
%           end
         VortPot_zm(k) = mean(VortPot(cc,k));

    end
    
           
  %% matrice [Nraggi*Ncerchi, totaltime]
  dVz_r_time (:,t) = dVz_dr(:);
  dVr_phi_time (:,t) = dVr_dphi(:);
  VortRelat_time(:,t) = VortRelat(:);
  VortPot_time(:,t) = VortPot(:);
  VortPot_zm_time(:,t) = VortPot_zm(:);

  
  % derivata seconda  
  %dVz_2_dr = dVz_2./(dr^2);
  
  %% matrice [Nraggi*Ncerchi, totaltime]
 % dVz_2_r_time (:,t) = dVz_2_dr(:);

end

% % %% -----------------------------------------------------------------------
% % salvataggio 
% %% ------------------------------------------------------------------------
% 
% fileout2 = ['dVz_dr_time_07'];
fileout2 = [roots,names{cexp,2},'/',output,'/dVz_dr_time_',num2str(Nraggi),'_',num2str(Ncerchi),'_',num2str(nFrames)];
filename = sprintf('%s.mtx',fileout2);
fid = fopen(filename,'wb');
fwrite(fid,size(dVz_r_time,1),'ulong');
fwrite(fid,size(dVz_r_time,2),'ulong');
fwrite(fid,dVz_r_time(:),'float');
fclose(fid);

% fileout3 = ['dVr_phi_time_07'];
fileout3 = [roots,names{cexp,2},'/',output,'/dVr_dphi_time_',num2str(Nraggi),'_',num2str(Ncerchi),'_',num2str(nFrames)];
filename = sprintf('%s.mtx',fileout3);
fid = fopen(filename,'wb');
fwrite(fid,size(dVr_phi_time,1),'ulong');
fwrite(fid,size(dVr_phi_time,2),'ulong');
fwrite(fid, dVr_phi_time(:),'float');
fclose(fid);

% fileout4 = ['VortPot_time_07'];
fileout4 = [roots,names{cexp,2},'/',output,'/VortPot_time_',num2str(Nraggi),'_',num2str(Ncerchi),'_',num2str(nFrames)];
filename = sprintf('%s.mtx',fileout4);
fid = fopen(filename,'wb');
fwrite(fid,size(VortPot_time,1),'ulong');
fwrite(fid,size(VortPot_time,2),'ulong');
fwrite(fid,VortPot_time(:),'float');
fclose(fid);

% 
% fileout5 = ['VortRelat_time_07'];
fileout5 = [roots,names{cexp,2},'/',output,'/VortRelat_time_',num2str(Nraggi),'_',num2str(Ncerchi),'_',num2str(nFrames)];
filename = sprintf('%s.mtx',fileout5);
fid = fopen(filename,'wb');
fwrite(fid,size(VortRelat_time,1),'ulong');
fwrite(fid,size(VortRelat_time,2),'ulong');
fwrite(fid,VortRelat_time(:),'float');
fclose(fid);

% fileout6 = ['VortPot_zm_time_07'];
fileout6 = [roots,names{cexp,2},'/',output,'/VortPot_zm_time_',num2str(Nraggi),'_',num2str(Ncerchi),'_',num2str(nFrames)];
filename = sprintf('%s.mtx',fileout6);
fid = fopen(filename,'wb');
fwrite(fid,size(VortPot_zm_time,1),'ulong');
fwrite(fid,size(VortPot_zm_time,2),'ulong');
fwrite(fid,VortPot_zm_time(:),'float');
fclose(fid);
%VortPot_zm_time = zeros (Ncerchi,totaltime);
% % filename = sprintf('%s.mtx',fileout2);
% % fid = fopen(filename,'wb');
% % fwrite(fid,size(dVz_2_r_time,1),'ulong');
% % fwrite(fid,size(dVz_2_r_time,2),'ulong');
% % fwrite(fid,dVz_2_r_time(:),'float');
% % fclose(fid);


