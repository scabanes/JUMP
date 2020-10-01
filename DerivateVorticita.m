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
% ##################################################################################################################################################
%                                                                                                                                              LOAD:
% ##################################################################################################################################################
[Vz_int,b] = loadmtx([roots,Name,NameVt]);
[Vr_int,b] = loadmtx([roots,Name,NameVr]);
%%----------------- mise en cm
Vz_int = Vz_int.*Tocm;
Vr_int = Vr_int.*Tocm;
%%----------------- Matrice initiated to zero
dVz_r_time = zeros (Nti*Nri,nTime);
dVr_r_time = zeros (Nti*Nri,nTime);
VortRelat_time = zeros (Nti*Nri,nTime);
VortPot_time = zeros (Nti*Nri,nTime);
VortPot_zm_time = zeros (Nri,nTime);
%%----------------- Time counter
iit=0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              TIME LOOP:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=Itime:Tmax
    clear Vz
    clear Vr
    disp (t);
    Vz = reshape (Vz_int(:,t),Nti,Nri);
    Vr = reshape (Vr_int(:,t),Nti,Nri);
        
  dVz = zeros(Nti,Nri);
  dVr = zeros(Nti,Nri);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VETTORE-matrice DVz/Dr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% ------------- Derivata sul raggio
  % calcolo le differenze sulla V zonale
    for i=2:Nri-1
    
       %%% dal punto r(2) uso le differenze finite centrate fino al
       %%% punto r(Ncerchi-1); per il punto r(1) dVz =0;
       % derivata prima
       dVz(:,i) = r(i+1)*Vz(:,i+1) - r(i-1)*Vz(:,i-1);
       
       % derivata seconda
     %   dVz_2(:,i) = Vz(:,i+1) - 2*Vz(:,i) + Vz(:,i-1);
    
    end
% derivata prima
  dVz_dr = dVz./(2*dr);
%--------------------------------------------------------------------------
% ------------- Derivata su theta
     % calcolo le differenze sulla V radiale
    for ii=2:Nti-1
    
       %%% dal punto r(2) uso le differenze finite centrate fino al
       %%% punto r(Ncerchi-1); per il punto r(1) dVz =0;
       
       % derivata prima
       dVr(ii,:) = Vr(ii+1,:) - Vr(ii-1,:);
       
       % derivata seconda
     %   dVz_2(:,i) = Vz(:,i+1) - 2*Vz(:,i) + Vz(:,i-1);
    
     end
% derivata prima
  dVr_dphi = dVr./(2*dtheta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRIMI e ULTIMI PUNTI SPAZIALE per derivate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% ------------- Derivata sul raggio
  for j=1:Nti
  aa = ~isnan(Vz(j,:));
  aa_min = min(find(aa==1));
  aa_max = max(find(aa==1));
  %find(aa==1)
    %%per il punto r(Ncerchi) uso le differenze finite all'indietro
  dVz_dr (j,aa_max) = (r(aa_max).*Vz(j,aa_max) - r(aa_max-1).*Vz(j,aa_max-1))./(dr);
    %%per il punto r(0) uso le differenze finite in avanti
  dVz_dr (j,aa_min) = (r(aa_min+1).*Vz(j,aa_min+1) - r(aa_min).*Vz(j,aa_min))./(dr);
  end
%--------------------------------------------------------------------------
% ------------- Derivata su theta
    for jj=1:Nri
  aaa = ~isnan(Vr(:,jj));
  aaa_min = min(find(aaa==1));
  aaa_max = max(find(aaa==1));
  %find(aa==1)
    %%per il punto r(Ncerchi) uso le differenze finite all'indietro
  dVr_dphi (aaa_max,jj) = (Vr(aaa_max,jj) - Vr(aaa_max-1,jj))./(dtheta);
    %%per il punto r(0) uso le differenze finite in avanti
  dVr_dphi (aaa_min,jj) = (Vr(aaa_min+1,jj) - Vr(aaa_min,jj))./(dtheta);
    end
%--------------------------------------------------------------------------
% ------------- Vorticita
    VortRelat = zeros(Nti,Nri);
    VortPot = zeros(Nti,Nri);
    VortPot1 = zeros(Nti,Nri);
    VortPot_zm = zeros(Nri,1);
    for k=1:Nri
         VortRelat(:,k) = (1./r(k))*(dVz_dr(:,k) - dVr_dphi(:,k));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% matrice [Nraggi*Ncerchi, totaltime]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% ------------- Condizione su r=0
% !!! This is important for not having Inf value.
% se r=0 exists then we put a nan instead of having infinity in the
% vorticity which is: vort ~ 1/r (...)
VortRelat(:,find(r==0))=nan;
VortPot(:,find(r==0))=nan;
VortPot_zm(find(r==0))=nan;
%--------------------------------------------------------------------------
% ------------- Save Matrices in time
iit=iit+1;
dVz_r_time (:,iit) = dVz_dr(:);
dVr_phi_time (:,iit) = dVr_dphi(:);
VortRelat_time(:,iit) = VortRelat(:);
VortPot_time(:,iit) = VortPot(:);
VortPot_zm_time(:,iit) = VortPot_zm(:);

  
  % derivata seconda  
  %dVz_2_dr = dVz_2./(dr^2);
  
  %% matrice [Nraggi*Ncerchi, totaltime]
 % dVz_2_r_time (:,t) = dVz_2_dr(:);

end
%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                  PLOTS
%  ###########################################################################################################################
%  ###########################################################################################################################


%  ###########################################################################################################################
%  ###########################################################################################################################

% % %% -----------------------------------------------------------------------
% % salvataggio 
% %% ------------------------------------------------------------------------
% 
% fileout2 = ['dVz_dr_time_07'];
fileout2 = [roots,Name,'/dVz_dr_time_',num2str(Nti),'_',num2str(Nri),'_',num2str(nTime)];
filename = sprintf('%s.mtx',fileout2);
fid = fopen(filename,'wb');
fwrite(fid,size(dVz_r_time,1),'ulong');
fwrite(fid,size(dVz_r_time,2),'ulong');
fwrite(fid,dVz_r_time(:),'float');
fclose(fid);

% fileout3 = ['dVr_phi_time_07'];
fileout3 = [roots,Name,'/dVr_dphi_time_',num2str(Nti),'_',num2str(Nri),'_',num2str(nTime)];
filename = sprintf('%s.mtx',fileout3);
fid = fopen(filename,'wb');
fwrite(fid,size(dVr_phi_time,1),'ulong');
fwrite(fid,size(dVr_phi_time,2),'ulong');
fwrite(fid, dVr_phi_time(:),'float');
fclose(fid);

% fileout4 = ['VortPot_time_07'];
fileout4 = [roots,Name,'/VortPot_time_',num2str(Nti),'_',num2str(Nri),'_',num2str(nTime)];
filename = sprintf('%s.mtx',fileout4);
fid = fopen(filename,'wb');
fwrite(fid,size(VortPot_time,1),'ulong');
fwrite(fid,size(VortPot_time,2),'ulong');
fwrite(fid,VortPot_time(:),'float');
fclose(fid);

% 
% fileout5 = ['VortRelat_time_07'];
fileout5 = [roots,Name,'/VortRelat_time_',num2str(Nti),'_',num2str(Nri),'_',num2str(nTime)];
filename = sprintf('%s.mtx',fileout5);
fid = fopen(filename,'wb');
fwrite(fid,size(VortRelat_time,1),'ulong');
fwrite(fid,size(VortRelat_time,2),'ulong');
fwrite(fid,VortRelat_time(:),'float');
fclose(fid);

% fileout6 = ['VortPot_zm_time_07'];
fileout6 = [roots,Name,'/VortPot_zm_time_',num2str(Nti),'_',num2str(Nri),'_',num2str(nTime)];
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


