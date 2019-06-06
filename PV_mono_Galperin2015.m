%function [LL_m,Mono_rp,Ndata] = jrh_L_mono_Torino(X,resolution);
%% Parametri da scegliere: -
clear all
close all
%**********************************************************
%           PARAMETRI DA SCEGLIERE
%**********************************************************
roots = '/media/simon/simon/ESP_29/'; % Root path..
run([roots,'InfosFile.m'])
%**********************************************************
Create_Grid_pol_Galperin2015
dr=R/Nri;
r=(1:Nri)*dr;
% ##################################################################################################################################################
%                                                                                                                                              LOAD:
% ##################################################################################################################################################
% Carichiamo input--------------------------
nFrames = Tmax-Itime+1;
[Pvpolrm,b]=loadmtx([roots,'/VortPot_time_',num2str(Nti),'_',num2str(Nri),'_',num2str(nFrames)]);
[Vz,b] = loadmtx([roots,NameVt]);

% b(2)=nFrames
% [Pvpolrm,b]=loadmtx('../VortPot_time_07');
% [Vz,b]=loadmtx('../Vz07_540_160_1419');
% [Pvpolrm,b]=loadmtx('../VortPot_540_160_2160');
% [Vz,b]=loadmtx('../Vz_540_160_2160');
% Create_Grid_pol_Torino;
% b(2)=5
% % % % criterio = 1; 
% % % pcolor(Grid_Xp_m, Grid_Yp_m,reshape(Pvpolrm(:,1),Nraggi,Ncerchi));shading interp
% % % title('Selectionner Rmin & Rmax')
% % % [Rmin,Rmax] = ginput(2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Description
%jrh_L_mono - Monotonizes column vector of evenly spaced observations, X,
%and outputs the L_m (analogous to Thorpe's Scale L_T [the RMS of the
%displacement for the displaced particles]), and the monotonized values in
%a column vector. If X is an array, this function will output the L_m for
%each column into a row vector and output the monotonized columns into an
%array the same size as X.
%
%
%Syntax: [L_m,Mono] = jrh_L_mono(X,scale,dir)
%
%Inputs:  X     - data to be monotonized
%         scale - distance between each observation (default is unit
%                 spacing)
%         dir   - profiles direction, mostly ascending or mostly
%                 descending. This should be the same for all profiles
%                 (default is mostly ascending)
%Outputs: L_m   - The monotonization (Thorpe) scale; zero is returned if no
%                 the values are already monotonized
%         Mono  - The monotonized values of the 
%
%Example:
%
%Created and tested with Matlab R2012a, February 2014

%Author: Jesse Hoemann, College of Marine Science, University of South
%        Florida
%email: jesseh1@mail.usf.edu
% Last revision: N/A (for now)
%
%Future revisions may account for cells of arrays to compare different data
%sets of varying size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ##################################################################################################################################################
%                                                                                                                                   NUMERO DEI DATI:
% ##################################################################################################################################################
Ndata_raggio_t = zeros(Nti,nFrames);
Ndata_cerchi_t = zeros(Nri,nFrames);
for it=1:nFrames
    pvPolr=reshape(Pvpolrm(:,it),Nti,Nri);

    for iphi=1:Nti
         % Salviamo il numero di dati su ogni raggi
         aa = ~isnan(pvPolr(iphi,:));
         Ndata_raggio(iphi) = length(pvPolr(iphi,aa));
    end
    for iR = 1:Nri
        % Salviamo il numero di dati su ogni cerchio
        kkk = ~isnan(pvPolr(:,iR));
        Ndata_cerchi(iR) = length(pvPolr(kkk,iR));
    end

% Numero dei dati  su raggi e cerchi.
Ndata_raggio_t(:,it) = Ndata_raggio;
Ndata_cerchi_t(:,it) = Ndata_cerchi;
end
figure
subplot(2,1,1);
plot(Ndata_raggio_t,'.')
ylabel('Number of data points in azimuth')
subplot(2,1,2);
plot(Ndata_cerchi_t,'.')
ylabel('Number of data points on the radius')
PtDaTogliereR = [];% one has to put '[]' if no point to take out 
PtDaTogliereT = [1:Ntmin,Ntmax:360];% one has to put '[]' if no point to take out 
%##########################################################################
%Create_Grid_pol_Torino_m
figure
pcolor(Grid_Xp_cm,Grid_Yp_cm,reshape(Pvpolrm(:,1),Nti,Nri)); shading interp

% ##################################################################################################################################################
%                                                                                                                                   LOOP SUL TEMPO:
% ##################################################################################################################################################
Mono_rpt = zeros(Nti,Nri,nFrames);
Lm_t = zeros(Nti,nFrames);
Mono_rpt_zm = zeros(Nri,nFrames);
%
VortPot_zm = nan(Nri,1);
VortPot_zmt = zeros(Nri,nFrames);
%
Vz_zm = nan(Nri,1);
Vz_zmt = zeros(Nri,nFrames);
%
dVortPot_zm_dr = nan(Nri,1);
dVortPot_zmt_dr= zeros(Nri,nFrames);

for it=1:nFrames
it
pvPolr=reshape(Pvpolrm(:,it),Nti,Nri);
Vz_t=reshape(Vz(:,it),Nti,Nri);
% Togliamo i punti con pocchi dati schelto in PtDaTogliere;
pvPolr(:,PtDaTogliereR) = nan;
Vz_t(:,PtDaTogliereR) = nan;
pvPolr(PtDaTogliereT,:) = nan;
Vz_t(PtDaTogliereT,:) = nan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=size(pvPolr);
LL_m = zeros(a(1),1);
Mono_rp = zeros(a);
    for iphi=1:Nti
        aa = ~isnan(pvPolr(iphi,:));
        pvPolr_NoNan = pvPolr(iphi,aa)';
    % How many data we have on the radius
    % it allows to define the angolar area in 
    % which we perform the zonal mean of the 
    % the monotonized profils.
    
        unsortI = reshape(1:length(pvPolr_NoNan),size(pvPolr_NoNan));%,size(X_NoNan);
    %[Mono,sortI] = sort(X_NoNan,1,'ascend');
        [Mono,sortI] = sort(pvPolr_NoNan,1,'descend');

    
        L = (unsortI-sortI).*dr;
        L_m = squeeze((sum(L.^2)./sum((L~=0))).^0.5);
   
        ii=1;
        jj=1;
        for iR = 1:Nri
            if(aa(jj) == 1)
                Mono_rp(iphi,iR) = Mono(ii);
                Mono_rp(iphi,iR) = Mono(ii);
                ii = ii+1;
            else
                Mono_rp(iphi,iR) = nan;
                Mono_rp(iphi,iR) = nan;
            end
            jj = jj+1;
        end


    LL_m(iphi) = L_m;
    clear L_M
    
    end

    for iR = 1:Nri
        % zonal mean della PV
        cc = ~isnan(pvPolr(:,iR));
        VortPot_zm(iR) = mean(pvPolr(cc,iR));
        % zonal mean della velocita zonale
        ccc = ~isnan(Vz_t(:,iR));
        Vz_zm(iR) = mean(Vz_t(ccc,iR));
    end
    
    pprima=1;
    for iR = 1:Nri
        if (~isnan(VortPot_zm(iR)) == 1 & pprima == 1) % derivata in avanti
            pprima =2;
            dVortPot_zm_dr(iR) = (VortPot_zm(iR+1)-VortPot_zm(iR))/dr;
        elseif (~isnan(VortPot_zm(iR)) == 1 & iR == Nri)% derivata all indietro
            dVortPot_zm_dr(iR) = (VortPot_zm(iR)-VortPot_zm(iR-1))/(dr);
        elseif (~isnan(VortPot_zm(iR)) == 1 & pprima == 2)% derivata centrata
            dVortPot_zm_dr(iR) = (VortPot_zm(iR+1)-VortPot_zm(iR-1))/(2.*dr);
        else
            dVortPot_zm_dr(iR) = nan;
        end
    end

% Facciamo una media su un settore angolare che scegliamo
% con indici presi su Ndata con un criterio arbitrario
% In realta non possiamo fare una media zonale sulla monotonizzazione
% i raggi non sono piu raggi dopo di avere monotonizzato. Quindi
% i dati della vorticita devono essere sugli stessi raggi
% per essere paragonati e mediati. Se il criterio e piccolo
% significa che ci son piu o meno lo stesso numero di dati 
% e una media non e troppo sbagliata.
% % % % idxData = find(Ndata(:) >= (max(Ndata(:))-criterio));
% % % % Mono_rp_zm = NaN(Ncerchi,1);
% % % % idxRmin = round(Rmin(2)/dr); 
% % % % idxRmax = round(Rmax(2)/dr); 
% % % % for iR = idxRmin:idxRmax%Ncerchi
% % % %     Mono_rp_pm = Mono_rp(idxData,iR);
% % % %     ff = ~isnan(Mono_rp_pm);
% % % %     Mono_rp_pm_SNan = Mono_rp_pm(ff);
% % % %     Mono_rp_zm(iR) = mean(Mono_rp_pm_SNan);
% % % % end


Lm_t(:,it) = LL_m;
Mono_rpt(:,:,it) = Mono_rp;

% PV & Velocita mediate
VortPot_zmt(:,it) = VortPot_zm;
Vz_zmt(:,it) = Vz_zm;
% Derivata della vorticita e media zonale
dVortPot_zmt_dr(:,it) = dVortPot_zm_dr;

% % % % Mono_rpt_zm(:,it) = Mono_rp_zm;
end
figure
pcolor(Grid_Xp_cm,Grid_Yp_cm,pvPolr); shading interp
% ##################################################################################################################################################
%                                                                                                                              VORTICITA POTENTIALE:
% ##################################################################################################################################################
% #########################################################################
% Facciamo diverse proceduri di medie temporalle et spaziale:
% Procedure I: medie temporale e spaziale sulla scala di Thorpe
%              La scala e fatta sui profili instatanei e dopo facciamo le
%              medie temporale e sul settore angolare scelto.
%              !!! c'e un criterio sulla griglia!!! --> The criterion on
%              the grid might be a bit subjective and make that procedure
%              less convicing.
% Procedure II: medie temporale e spaziale sui profili di vorticit�
%               Facciamo la media zonale e temporale dei profili di
%               vorticit�. Dopo monotoniziamo e esce la scala su un solo 
%               profilo mediato. --> Averaging procedure are likely to kill
%               turbulence properties
% Procedure III: media spaziale sulla vorticit� e temporale sulla scala 
%               Facciamo la media zonale della vorticit�, monotoniziamo
%               e esce la scala di Thorpe per ogni instantaneo. Dopo
%               mediamo la scala sul tempo. ----> THE ONE I WOULD PREFERE
% #########################################################################
% #########################################################################
%                       Medie temporale & Spatiale
% #########################################################################
%%                             Procedure I
% #########################################################################
% criterio che taglia i profili radiali troppo corti per avere la
% risoluzione della scala di Thorpe, i.e. ai bordi angolari della
% camera.
criterio = 10;
disp('Attenzione la scala di Thorpe deve essere al massimo la met� di:')
disp(num2str(dr*criterio))
% Media temporale e seleziona un settore definito 
% dove ci sono meno dati del criterio
Ndata_mt=mean(Ndata_raggio_t');
for idxphi=1:Nti
%    if(idxphi >= idxCmin & idxphi <= idxCmax)
    if(Ndata_mt(idxphi) >= criterio)
        sss=~isnan(Lm_t(idxphi,:));
        Lm_mt(idxphi) = mean(Lm_t(idxphi,sss)');
    else
        Lm_mt(idxphi) = NaN;
    end
end
% Media sull settore
ddd=~isnan(Lm_mt);
Lm_mtzm = mean(Lm_mt(ddd));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        THORPE SCALE
disp('scala di Thorpe I:')
disp(num2str(Lm_mtzm))
% #########################################################################
%%                             Procedure II
% #########################################################################
VortPot_zmtm = mean(VortPot_zmt');
ggg = ~isnan(VortPot_zmtm);
unsortII = [1:length(VortPot_zmtm(ggg))];%,size(X_NoNan);
[MonoII,sortII] = sort(VortPot_zmtm(ggg)',1,'descend');
LII = (unsortII'-sortII).*dr;
LII_m = squeeze((sum(LII.^2)./sum((LII~=0))).^0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        THORPE SCALE
disp('scala di Thorpe II:')
disp(num2str(LII_m))

figure; hold on
plot(VortPot_zmtm,r,'.b','Linewidth',2)
plot(MonoII,r(ggg),'k','Linewidth',1.5)
xlabel('Potential Vorticity (m^-1 s^-1)')
ylabel('Radius (m)')
title(['Zonal and temporal average of the PV --> PII'])%--> L_M = ',num2str(LII_m),' (m).'])
legend('Potential vorticity','Monotonized vorticity')

% #########################################################################
%%                            Procedure  III
% #########################################################################
for it=1:nFrames
VortPot_zm1t=VortPot_zmt(:,it)';
fff = ~isnan(VortPot_zm1t);
unsortIII = [1:length(VortPot_zm1t(fff))];%,size(X_NoNan);
[MonoIII,sortIII] = sort(VortPot_zm1t(fff)',1,'descend');
LIII = (unsortIII'-sortIII).*dr;
LIII_m = squeeze((sum(LIII.^2)./sum((LIII~=0))).^0.5);
LIIIm_t(it)= LIII_m;
end
LIIIm_mt=mean(LIIIm_t);

figure; hold on
plot(VortPot_zm1t,r,'.b','Linewidth',2)
plot(MonoIII,r(fff),'k','Linewidth',1.5)
xlabel('Potential Vorticity (m^{-1}.s^{-1})')
ylabel('Radius (m)')
title(['Instantaneous profiles of the PV zonal mean --> PIII'])% & <L_M>_T = ',num2str(LIIIm_mt),' (m).'])
legend('Potential vorticity','Monotonized vorticity')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        THORPE SCALE
disp('scala di Thorpe III:')
disp(num2str(LIIIm_mt))

% ##################################################################################################################################################
%                                                                                               VORTICITA POTENTIALE + DERIVATA DELLA PV + VELOCITA:
% ##################################################################################################################################################
% ------------------------------- Fiure of profiles 1
figure; hold on
line(VortPot_zm1t,r,'Color','k','Linewidth',1.5)
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';
xlabel('PV m^{-1}.s^{-1}')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(Vz_zmt(:,it),r,'Parent',ax2,'Color','r','Linewidth',1.5)
line([0 0],get(ax2,'YLim'),'Parent',ax2,'Color','r')
ax2.XColor = 'r';
ax2.YColor = 'r';
xlabel('zonal velocity m.s^{-1}')
% ------------------------------- Fiure of profiles 2
figure; hold on
line(VortPot_zm1t,r,'Color','k','Linewidth',1.5)
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';
xlabel('PV m^{-1}.s^{-1}')
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
line(dVortPot_zmt_dr(:,it),r,'Parent',ax2,'Color','r','Linewidth',1.5)
line([0 0],get(ax2,'YLim'),'Parent',ax2,'Color','r')
ax2.XColor = 'r';
ax2.YColor = 'r';
xlabel('DPV/Dr m^{-1}.s^{-1}')




% ##################################################################################################################################################
%                                                                                                                                          SALVIAMO:
% ##################################################################################################################################################
% jrh_savemtx(['Lms'],Lm_t);
% jrh_savemtx(['Ndata_t'],Ndata_t);
% jrh_savemtx(['pvMonoPolrs'],reshape(Mono_rpt,Nraggi*Ncerchi,[]));
% jrh_savemtx(['Mono_rpt_zm'],Mono_rpt_zm);


