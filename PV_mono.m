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
%                                                                                                                                              LOAD:
% ##################################################################################################################################################
% Carichiamo input--------------------------
[Pvpolrm,b]=loadmtx([roots,Name,'/VortPot_time_',num2str(Nti),'_',num2str(Nri),'_',num2str(nTime)]);
% [PvRelat,b]=loadmtx([roots,Name,'/VortRelat_time_',num2str(Nti),'_',num2str(Nri),'_',num2str(nTime)]);
[Vz,b] = loadmtx([roots,Name,NameVt]);
% ##################################################################################################################################################
%                                                                                                                                   NUMERO DEI DATI:
% ##################################################################################################################################################
Ndata_raggio_t = zeros(Nti,nTime);
Ndata_cerchi_t = zeros(Nri,nTime);
for it=1:nTime
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
if(0==1)
figure
subplot(2,1,1);
plot(Ndata_raggio_t,'.')
ylabel('Number of data points in azimuth')
subplot(2,1,2);
plot(Ndata_cerchi_t,'.')
ylabel('Number of data points on the radius')
end
%##########################################################################
figure
pcolor(Grid_Xp_cm,Grid_Yp_cm,reshape(Pvpolrm(:,1),Nti,Nri)); shading interp

%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   PV MONOTISAZIONE
%  ###########################################################################################################################
%  ###########################################################################################################################
Mono_rpt = zeros(Nti,Nri,nTime);
Lm_t = zeros(Nti,nTime);
Mono_rpt_zm = zeros(Nri,nTime);
%
VortPot_zm = nan(Nri,1);
VortPot_zmt = zeros(Nri,nTime);
%
Vz_zm = nan(Nri,1);
Vz_zmt = zeros(Nri,nTime);
%
dVortPot_zm_dr = nan(Nri,1);
dVortPot_zmt_dr= zeros(Nri,nTime);

% ##################################################################################################################################################
%                                                                                                                                   LOOP SUL TEMPO:
% ##################################################################################################################################################
iit=0;
for it=Itime:Tmax
disp(it)
%--------------------------------------------------------------------------
% ------------- Reshape
pvPolr=reshape(Pvpolrm(:,it),Nti,Nri);
Vz_t=reshape(Vz(:,it),Nti,Nri);
%--------------------------------------------------------------------------
% ------------- Taglio sui punti
% Togliamo i punti con pocchi dati schelto in PtDaTogliere;
pvPolr(:,PtDaTogliereR) = nan;
Vz_t(:,PtDaTogliereR) = nan;
pvPolr(PtDaTogliereT,:) = nan;
Vz_t(PtDaTogliereT,:) = nan;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                INSTATANEOUS OF PV MONO:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=size(pvPolr);
LL_m = zeros(a(1),1);
Mono_rp = zeros(a);
    for iphi=1:Nti%----------------------------------------------------------------------Loop on theta
        %--------------------------------------------------------------------------
        % ------------- Sorting instantaneoous of PV with nonan.
        aa = ~isnan(pvPolr(iphi,:));
        pvPolr_NoNan = pvPolr(iphi,aa)';
    
        unsortI = reshape(1:length(pvPolr_NoNan),size(pvPolr_NoNan));%,size(X_NoNan);
        %[Mono,sortI] = sort(X_NoNan,1,'ascend');
        [Mono,sortI] = sort(pvPolr_NoNan,1,'descend');
 
        L = (unsortI-sortI).*dr;
        L_m = squeeze((sum(L.^2)./sum((L~=0))).^0.5);
        %--------------------------------------------------------------------------
        % ------------- Save Monotonized prf in a Martix. Non-used.
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
        %--------------------------------------------------------------------------
        % ------------- Save LM for PI

        LL_m(iphi) = L_m;
        clear L_M
    end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                   ZONAL MEAN ON THE PV:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iR = 1:Nri%----------------------------------------------------------------------Loop on radius
        % zonal mean della PV
        cc = ~isnan(pvPolr(:,iR));
        VortPot_zm(iR) = mean(pvPolr(cc,iR));
        % zonal mean della velocita zonale
        ccc = ~isnan(Vz_t(:,iR));
        Vz_zm(iR) = mean(Vz_t(ccc,iR));
    end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        DERIVATA DELLA PV IN ZONAL MEAN:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
iit=iit+1;
%--------------------------------------------------------------------------
% ------------- For PI
Lm_t(:,iit) = LL_m;
Mono_rpt(:,:,iit) = Mono_rp;
%--------------------------------------------------------------------------
% ------------- For PII & PIII
% PV & Velocita mediate
VortPot_zmt(:,iit) = VortPot_zm;
%--------------------------------------------------------------------------
% ------------- AUTRES
Vz_zmt(:,iit) = Vz_zm;
% Derivata della vorticita e media zonale
dVortPot_zmt_dr(:,iit) = dVortPot_zm_dr;

% % % % Mono_rpt_zm(:,it) = Mono_rp_zm;
end
figure
pcolor(Grid_Xp_cm,Grid_Yp_cm,pvPolr); shading interp

%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   PROCEDURE
%  ###########################################################################################################################
%  ###########################################################################################################################

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
% Procedure II: medie temporale e spaziale sui profili di vorticita
%               Facciamo la media zonale e temporale dei profili di
%               vorticita. Dopo monotoniziamo e esce la scala su un solo 
%               profilo mediato. --> Averaging procedure are likely to kill
%               turbulence properties
% Procedure III: media spaziale sulla vorticita e temporale sulla scala 
%               Facciamo la media zonale della vorticita, monotoniziamo
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
criterio = 30;
disp('Attenzione la scala di Thorpe deve essere al massimo la meta di:')
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
disp('scala di Thorpe I - Full instantaneous (used):')
disp(num2str(Lm_mtzm))


%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                  PLOTS
%  ###########################################################################################################################
%  ###########################################################################################################################
% #################################################################################################################
%                                                                                                   Letter's Plots:
% #################################################################################################################
% We set plots to write a letter
iitt = 81; % 90 % choose a time to plot here
Hmean = 4 ; % hauteur moyenne de fluide, 4 cm.
% #########################################################################
%                                                                       PV:
% #########################################################################
%__________________________________________________________________________
%                           Plots Config
%__________________________________________________________________________
figure
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
   'Color',[1 1 1],'PaperPositionMode','auto')

axes('FontSize',18,'Linewidth',2,'FontName','times',...
     'TickLength',[0.01; 0.03],...%'Xdir','reverse',...
     'Position',[0.13 0.17 0.82 0.75])   
hold on
subplot ('Position',[0.13 0.17 0.82 0.75]) 
hold on
%__________________________________________________________________________
%----------------------------------Plot QGPV

%-----------------------PLOTS instantaneous non-averaged data
pvPolr=reshape(Pvpolrm(:,iitt),Nti,Nri);
plot(pvPolr(260,:)/(Omega/Hmean),r,'-k','linewidth',2.)
plot(Mono_rpt(260,:,iitt)/(Omega/Hmean),r,'-r','linewidth',1.)
box on
% ylim ([11 27])
ylim ([9 27])


% #########################################################################
%                                                                  U Zonal:
% #########################################################################
%__________________________________________________________________________
%                           Plots Config
%__________________________________________________________________________
figure
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
   'Color',[1 1 1],'PaperPositionMode','auto')

axes('FontSize',18,'Linewidth',2,'FontName','times',...
     'TickLength',[0.01; 0.03],...%'Xdir','reverse',...
     'Position',[0.13 0.17 0.82 0.75])   
hold on
subplot ('Position',[0.13 0.17 0.82 0.75]) 
hold on
%__________________________________________________________________________
%----------------------------------Plot QGPV

%-----------------------PLOTS instantaneous non-averaged data
plot(Vz_zmt(:,iitt),r,'-k','linewidth',2.)
box on
% ylim ([11 27])
ylim ([9 27])

% #################################################################################################################
%                                                                                                    Letter's Maps:
% #################################################################################################################

% #########################################################################
%                                                       Vorticita relativa:
% #########################################################################
% Path to colorbar repository
addpath('/home/simon/Bureau/Esperimento-DICEA/JUMP/Altre-funzioni')
addpath('/home/simon/Bureau/Esperimento-DICEA/JUMP/Altre-funzioni/cbrewer')
%
[PvRelat,b]=loadmtx([roots,Name,'/VortRelat_time_',num2str(Nti),'_',num2str(Nri),'_',num2str(nTime)]);
pvPolr=reshape(PvRelat(:,iitt),Nti,Nri);
fig=figure
pcolor(GridR(Ntmin:Ntmax,:).*cmpx,GridT(Ntmin:Ntmax,:),pvPolr(Ntmin:Ntmax,:)); shading interp
% contourf(GridR(180:end,:),GridT(180:end,:),pvPolr(180:end,:),50)
colormap(cbrewer('div', 'RdBu', 30))
colorbar
mycmap = get(fig,'Colormap')
set(fig,'Colormap',flipud(mycmap))
caxis([-3 3]);
title('Potential vorticity')

% #########################################################################
%                                                       Zonal and
%                                                       meridional
%                                                       velocity:
% #########################################################################
% --------- zonal velocity
[Vtheta_tot,b]=loadmtx([roots,Name,NameVt]);
Vtheta=reshape(Vtheta_tot(:,iitt),Nti,Nri);
fig=figure
pcolor(GridR(Ntmin:Ntmax,:).*cmpx,GridT(Ntmin:Ntmax,:),Vtheta(Ntmin:Ntmax,:)); shading interp
% contourf(GridR(180:end,:),GridT(180:end,:),pvPolr(180:end,:),50)
colormap(cbrewer('div', 'RdBu', 30))
colorbar
mycmap = get(fig,'Colormap')
set(fig,'Colormap',flipud(mycmap))
caxis([-3 3]);
title('Zonal velocity')

% --------- radial velocity
[Vr_tot,b]=loadmtx([roots,Name,NameVr]);
Vr=reshape(Vr_tot(:,iitt),Nti,Nri);
fig=figure
pcolor(GridR(Ntmin:Ntmax,:).*cmpx,GridT(Ntmin:Ntmax,:),Vr(Ntmin:Ntmax,:)); shading interp
% contourf(GridR(180:end,:),GridT(180:end,:),pvPolr(180:end,:),50)
colormap(cbrewer('div', 'RdBu', 30))
colorbar
mycmap = get(fig,'Colormap')
set(fig,'Colormap',flipud(mycmap))
caxis([-1 1]);
title('radial velocity')
% #################################################################################################################
%                                                                                                       Multiplots:
% #################################################################################################################
cmtom = 100; % convert from cm to m
%__________________________________________________________________________
%                           Plots Config
%__________________________________________________________________________
figure; hold on
ax1 = axes('FontSize',18,'Linewidth',2,'FontName','times',...
     'TickLength',[0.01; 0.03])   
%_________________________________________________________________________
%----------------------------------Plot QGPV
pvPolr=reshape(Pvpolrm(:,iitt),Nti,Nri);
line(pvPolr(260,:).*Hmean,r,'Color','k','Linewidth',2.5)
%----------------------------------Plot sorted QGPV
line(Mono_rpt(260,:,iitt).*Hmean,r,'Color','r','Linewidth',1.5)
% ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';
% xlabel('PV m^{-1}.s^{-1}')
ax1_pos = ax1.Position; % position of first axes
% ylim ([10 25])
% xlim ([5 25])
ylim ([10 23])
xlim ([7 25])
%_________________________________________________________________________
%----------------------------------Plot Zonal velocity
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none','FontSize',18,'Linewidth',2,'FontName','times',...
     'TickLength',[0.01; 0.03]);
line(Vz_zmt(:,iitt)./cmtom,r,'Parent',ax2,'Color','k','Linewidth',2.5)
ax2.XColor = 'k';
ax2.YColor = 'k';

% ylim ([10 25])
% xlim ([-2 3.6]./cmtom)
ylim ([10 23])
xlim ([-2 0.5]./cmtom)
title('Pv in s^{-1} and Vz en m s^{-1}')
%% ###########################################################################################################################
%  ###########################################################################################################################
%  ###########################################################################################################################
%  ###########################################################################################################################








if(1==0)
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
for it=1:nTime
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

end


% ##################################################################################################################################################
%                                                                                                                                          SALVIAMO:
% ##################################################################################################################################################
% jrh_savemtx(['Lms'],Lm_t);
% jrh_savemtx(['Ndata_t'],Ndata_t);
% jrh_savemtx(['pvMonoPolrs'],reshape(Mono_rpt,Nraggi*Ncerchi,[]));
% jrh_savemtx(['Mono_rpt_zm'],Mono_rpt_zm);


