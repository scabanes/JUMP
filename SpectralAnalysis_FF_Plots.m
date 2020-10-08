load([roots,Name,'/SpectralAnalysis_FF_infos.mat'])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              LOAD DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[E_ky_t,a] = loadmtx([roots,Name,'/E_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
[EZ_ky_t,a] = loadmtx([roots,Name,'/EZ_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
[ER_ky_t,a] = loadmtx([roots,Name,'/ER_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
[Ex_kx_t,a] = loadmtx([roots,Name,'/Ex_kx_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
[Ey_kx_t,a] = loadmtx([roots,Name,'/Ey_kx_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
[Ex_ky_t,a] = loadmtx([roots,Name,'/Ex_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
[Ey_ky_t,a] = loadmtx([roots,Name,'/Ey_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
% [Emn_t,a] = loadmtx([roots,Name,'/Emn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
% close all
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                            SPECTRAl parameters & Typical length scales:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------------------------------- Length of the domain
Lx = x(end)-x(1); % en m
Ly = y(end)-y(1); % en m
dy = mean(y(2:end)-y(1:end-1));
%----------------------------------------------------------- Wavenumbers
kx_m = (2*pi*kx(1:Nkx))./Lx; % en rad.m-1
ky_m = (2*pi*ky(1:Nky))./Ly; % en rad.m-1
% ky_m = (2.*pi.*ky(1:Nky))./(2.*Ly); % en m-1
% kx_m = (2.*pi.*kx(1:Nkx))./(2.*Lx); % en m-1

%----------------------------------------------------------- Typical scales.
% j ai confiance dans ces echelles et le facteur 2pi.
iTab=0;
EKE_S = sum(mean(ER_ky_t,2));
iTab = iTab+1;
k_beta(iTab) = ((Cz/Ck).^(3/10)).*((beta.^3)/epsilon).^(1./5.);
L_beta(iTab) = (2*pi)./k_beta(iTab);
k_Rhines(iTab) = (beta/(2.*sqrt(2.*EKE_S))).^(1./2.);
L_Rhines(iTab) = (2*pi)./k_Rhines(iTab);
k_woods(iTab) = (0.24./((24/55)*Ck)).^(3./4.) .* (((2.*(Omega)).^3)./epsilon).^(1./2.);
L_woods(iTab) = ((2*pi)./k_woods(iTab));
k_RO(iTab) = (((2.*(Omega)).^3)./epsilon).^(1./2.);
L_RO(iTab) = ((2.*pi)./k_RO(iTab));
R_beta(iTab) = L_Rhines(iTab)/L_beta(iTab);
RO(iTab) = sqrt(2.*TKE)./(2.*Omega*2.5);
disp(['Typical scales :'])
disp(['L_beta = ',num2str(L_beta(iTab)),' m'])
disp(['L_Rhines = ',num2str(L_Rhines(iTab)),' m'])
disp(['L_woods = ',num2str(L_woods(iTab)),' m'])
disp(['R_beta = ',num2str(R_beta(iTab))])
disp(['RO = ',num2str(RO(iTab))])
% L_mag = 1.2 % magnet size in cm
% L_jets = 7 % jets typical size in cm
% n_mag=R/L_mag;
% n_jets=R/L_jets;
%----------------------------------------------------------- End typical scales
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Total Enegy:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0==1)
% -------------------------------- plots spectra
fig=figure; 
hax=axes; 
%
loglog(ky_m,mean(E_ky_t,2),'k', 'Linewidth',3)
% ylim([10.^(-10.) 10.^(-5.)])
% xlim([0 Nky-1])
box on
xlabel('ky [cm^{-1}]','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E Total (cm^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')
set(gca,'FontSize',13,'FontWeight','bold')
% title(['\epsilon = ',num2str(epsilon),' cm^2.s^{-3}', 'and \Omega = ',num2str(Omega),' rad/s'],'FontSize',11)
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                            Zonostrophy:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0==0)
% les quantites theorique sont divisees par L pour passer des m^3 s-2 a
% des m^2 s-2
% -------------------------------- Theoretical spectra
EZ_Theo = (Cz*beta.^2.).*ky_m(2:end).^(-5.); % [E] = m^3 s-2
ER_Theo =  Ck.*(epsilon.^(2./3.)).*ky_m(2:end).^(-5./3.); % [E] = m^3 s-2
% -------------------------------- plots spectra zonostrophic formulation
fig=figure; 
hax=axes; 
% loglog([0:Nky-1],mean(EZ_ky_t,2),'r')
loglog(ky(1:Nky),mean(EZ_ky_t,2),'r','LineWidth',3)
hold on
% loglog([0:Nky-1],mean(ER_ky_t,2),'k','LineWidth',2)
loglog(ky(1:Nky),mean(ER_ky_t,2),'k','LineWidth',3)
% loglog([0:Nky-1],mean(ER_ky_t+EZ_ky_t,2),'k','LineWidth',2)
% hold on
% loglog([1:Nky-1],EZ_Theo/(2*Ly),'--r') % [E] = m^2 s-2
% loglog([1:Nky-1],ER_Theo/(2*Ly),'--k') % [E] = m^2 s-2
loglog(ky(2:Nky),EZ_Theo./Ly,'--r','LineWidth',2) % [E] = m^2 s-2
loglog(ky(2:Nky),ER_Theo./Ly,'--k','LineWidth',2) % [E] = m^2 s-2
xlabel ('$k$ [rad.m$^{-1}$]','FontSize',18,'FontName','times','Interpreter','Latex')
ylabel ('$E $(m$^{2}$s$^{-2}$)','FontSize',18,'FontName','times','Interpreter','Latex')
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')
set(gca,'FontSize',13,'FontWeight','bold')
title(['Zonostrophic quantities : $\varepsilon$ = ',num2str(epsilon),' m$^2$.s$^{-3}$', 'and $\Omega$ = ',num2str(Omega),' rad/s'],'FontSize',18,'FontName','times','Interpreter','Latex')
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                   QNSE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  QNSE - Transverse(ky):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0==0)
EQNSE_A = (24/55)*Ck.*(epsilon.^(2./3.)).*ky_m(2:end).^(-5./3.); % [E] = m^3 s-2 
EQNSE_B = 0.24.*((2.*Omega).^2.).*(ky_m(2:end)).^(-3.); % [E] = m^3 s-2
EQNSE = EQNSE_A + EQNSE_B;
% -------------------------------- plots spectra 
fig=figure; 
hax=axes; 
loglog(ky_m,mean(Ex_ky_t,2),'k', 'Linewidth',3)
hold on
loglog(ky_m(2:end),EQNSE_B./Ly,'-b') % [E] = m^2 s-2
loglog(ky_m(2:end),EQNSE_A./Ly,'-b') % [E] = m^2 s-2
% loglog(ky_m(2:end),EQNSE,'--','Color','r', 'Linewidth',3.5) % [E] = m^2 s-2
%
loglog(ky_m(2:end),EZ_Theo./Ly,'--r','LineWidth',2) % [E] = m^2 s-2
% ylim([5.*10.^(-11.) 10.^(-5.)])
% xlim([ky_m(1) ky_m(end)])
box on
xlabel ('$k$ [rad.m$^{-1}$]','FontSize',18,'FontName','times','Interpreter','Latex')
ylabel ('$E $(m$^{2}$s$^{-2}$)','FontSize',18,'FontName','times','Interpreter','Latex')
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')
set(gca,'FontSize',13,'FontWeight','bold')
title(['E$_{Transverse}$(ky) : $\varepsilon$ = ',num2str(epsilon),' m$^2$.s$^{-3}$', 'and $\Omega$ = ',num2str(Omega),' rad/s'],'FontSize',18,'FontName','times','Interpreter','Latex')
saveas(fig,[roots,Name,'/Figures/Spectrum_ETrans_ky.svg']);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  QNSE - Longitudinal(ky):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0==1)
EQNSE_A = (18/55)*Ck.*(epsilon.^(2./3.)).*ky_m(2:end).^(-5./3.); 
EQNSE_B = 0.0926.*((2.*Omega/(2*pi)).^2.).*(ky_m(2:end)).^(-3.);
EQNSE = EQNSE_A + EQNSE_B;
% -------------------------------- plots spectra 
fig=figure; 
hax=axes; 
%
loglog(ky_m,mean(Ey_ky_t,2),'k', 'Linewidth',3)
hold on
loglog(ky_m(2:end),EQNSE_B/Ly,'-b') % [E] = m^2 s-2
loglog(ky_m(2:end),EQNSE_A/Ly,'-b') % [E] = m^2 s-2
loglog(ky_m(2:end),EQNSE/Ly,'--','Color','r', 'Linewidth',3.5) % [E] = m^2 s-2
ylim([5.*10.^(-11.) 10.^(-5.)])
xlim([ky_m(1) ky_m(end)])
box on
xlabel('ky [m^{-1}]','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E_1 (m^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')
set(gca,'FontSize',13,'FontWeight','bold')
title(['E_{Longitudinal}(ky) : \epsilon = ',num2str(epsilon),' m^2.s^{-3}', 'and \Omega = ',num2str(Omega),' rad/s'],'FontSize',11)
saveas(fig,[roots,Name,'/Figures/Spectrum_ELong_ky.svg']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  QNSE - Tansvere(kx):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0==1)
EQNSE_A = (24/55)*Ck.*(epsilon.^(2./3.)).*kx_m(2:end).^(-5./3.); % [E] = m^3 s-2 
EQNSE_B = 0.24.*((2.*Omega/(2*pi)).^2.).*(kx_m(2:end)).^(-3.); % [E] = m^3 s-2
EQNSE = EQNSE_A + EQNSE_B;
% -------------------------------- plots spectra 
fig=figure; 
hax=axes; 
loglog(kx_m,mean(Ey_kx_t,2),'k', 'Linewidth',3)
hold on
loglog(kx_m(2:end),EQNSE_B/Lx,'-b') % [E] = m^2 s-2
loglog(kx_m(2:end),EQNSE_A/Lx,'-b') % [E] = m^2 s-2
loglog(kx_m(2:end),EQNSE/Lx,'--','Color','r', 'Linewidth',3.5) % [E] = m^2 s-2
ylim([5.*10.^(-11.) 10.^(-5.)])
xlim([kx_m(1) kx_m(end)])
box on
xlabel('kx [m^{-1}]','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E_1 (m^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')
set(gca,'FontSize',13,'FontWeight','bold')
title(['E_{Tansvere}(kx) : \epsilon = ',num2str(epsilon),' m^2.s^{-3}', 'and \Omega = ',num2str(Omega),' rad/s'],'FontSize',11)
saveas(fig,[roots,Name,'/Figures/Spectrum_ETrans_kx.svg']);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  QNSE - Longitudinal(kx):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0==1)
EQNSE_A = (18/55)*Ck.*(epsilon.^(2./3.)).*kx_m(2:end).^(-5./3.); 
EQNSE_B = 0.0926.*((2.*Omega/(2*pi)).^2.).*(kx_m(2:end)).^(-3.);
EQNSE = EQNSE_A + EQNSE_B;
% -------------------------------- plots spectra 
fig=figure; 
hax=axes; 
loglog(kx_m,mean(Ex_kx_t,2),'k', 'Linewidth',3)
hold on
loglog(kx_m(2:end),EQNSE_B/Lx,'-b') % [E] = m^2 s-2
loglog(kx_m(2:end),EQNSE_A/Lx,'-b') % [E] = m^2 s-2
loglog(kx_m(2:end),EQNSE/Lx,'--','Color','r', 'Linewidth',3.5) % [E] = m^2 s-2
ylim([5.*10.^(-11.) 10.^(-5.)])
xlim([kx_m(1) kx_m(end)])
box on
xlabel('kx [m^{-1}]','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E_1 (m^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')
set(gca,'FontSize',13,'FontWeight','bold')
title(['E_{Longitudinal}(kx) : \epsilon = ',num2str(epsilon),' m^2.s^{-3}', 'and \Omega = ',num2str(Omega),' rad/s'],'FontSize',11)
saveas(fig,[roots,Name,'/Figures/Spectrum_ELong_kx.svg']);
end
% 
% disp(['#--> Here comes some quantities of ',Name,' :'])
% disp(['epsilon = ',num2str(epsilon),' cm^2.s^{-3}'])
% disp(['L_beta = ',num2str(L_beta)])