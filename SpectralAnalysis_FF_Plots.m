load([roots,Name,'/SpectralAnalysis_FF_infos.mat'])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              LOAD DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[EZ_ky_t,a] = loadmtx([roots,Name,'/EZ_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
[ER_ky_t,a] = loadmtx([roots,Name,'/ER_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
[E1_kx_t,a] = loadmtx([roots,Name,'/E1_kx_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
[E1_ky_t,a] = loadmtx([roots,Name,'/E1_ky_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
% [Emn_t,a] = loadmtx([roots,Name,'/Emn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
% close all
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                SPECTRA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx = x(end)-x(1); % en m
Ly = y(end)-y(1); % en m
dy = mean(y(2:end)-y(1:end-1));
R=2.5;
ky = (2*pi*[0:Nky-1])/Ly; % en m-1
kx = (2*pi*[0:Nkx-1])/Lx; % en m-1

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Zonostrophy:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% les quantites theorique sont divisees par Ro pour passer des m^3 s-2 a
% des m^2 s-2
% -------------------------------- Theoretical spectra
EZ_Theo = (Cz*beta.^2.).*ky(2:end).^(-5.); % [E] = m^3 s-2
ER_Theo = Ck.*(epsilon.^(2./3.)).*ky(2:end).^(-5./3.); % [E] = m^3 s-2
% ---------------------------------scales
L_beta=(epsilon./beta.^3.).^(1./5.);
% Lo_beta = (Ck/Cz).^(3./10.).*(J_root(1,1:Nrk)./J_root(2,1:Nrk)).^(1./2.);
% L_R = sqrt(2.*Urms./beta); 
% n_r=R/L_R;
% L_mag = 1.2 % magnet size in cm
% L_jets = 7 % jets typical size in cm
% n_mag=R/L_mag;
% n_jets=R/L_jets;
% -------------------------------- plots spectra zonostrophic formulation
fig=figure; 
hax=axes; 
loglog([0:Nky-1],mean(EZ_ky_t,2),'r')
hold on
loglog([0:Nky-1],mean(ER_ky_t,2),'k','LineWidth',2)
% loglog([0:Nky-1],mean(ER_ky_t+EZ_ky_t,2),'k','LineWidth',2)
% hold on
% loglog([1:Nky-1],EQNSE_ky/Ro,'--r')
loglog([1:Nky-1],EZ_Theo/Ro,'--r') % [E] = m^2 s-2
loglog([1:Nky-1],ER_Theo/Ro,'--k') % [E] = m^2 s-2

close all
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  QNSE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EQNSE_1 = (18/55)*Ck.*(epsilon.^(2./3.)).*kx(2:end).^(-5./3.); 
% EQNSE_2 = 0.0926.*((2.*omega.*sin(pi/3)).^2.).*(kx(2:end)).^(-3.);
% sur ky
EQNSE_1 = 0.626*(epsilon.^(2./3.)).*ky(2:end).^(-5./3.); % [E] = m^3 s-2 
EQNSE_2 = 0.24.*((2.*omega/(2*pi)).^2.).*(ky(2:end)).^(-3.); % [E] = m^3 s-2
EQNSE = EQNSE_1 + EQNSE_2;
% -------------------------------- plots spectra zonostrophic formulation
fig=figure; 
hax=axes; 
% loglog([0:Nkx-1],mean(E1_kx_t,2),'k')
% hold on
% loglog([1:Nkx-1],EQNSE_2/R,'-b')
% loglog([1:Nkx-1],EQNSE_1/R,'-b')
%
loglog([0:Nky-1],mean(E1_ky_t,2),'k')
hold on
loglog([1:Nky-1],EQNSE_2/Ro,'-b') % [E] = m^2 s-2
loglog([1:Nky-1],EQNSE_1/Ro,'-b') % [E] = m^2 s-2
loglog([1:Nky-1],EQNSE/Ro,'-r') % [E] = m^2 s-2






% ylim([10.^(-7.) 0.3].*cm2tom2)
% xlim([n(1) n(end)])
% line([R/Lhat_beta(nhat) R/Lhat_beta(nhat)],get(hax,'YLim'),'Color',[0 0 0])
% line([n_r n_r],get(hax,'YLim'),'Color',[0 0 0])
% line([n_mag n_mag],get(hax,'YLim'),'Color',[0 0 0])
% line([2*n_mag 2*n_mag],get(hax,'YLim'),'Color',[0 0 0])

% box on
% xlabel('n','FontSize',13,'FontWeight','bold','Color','k')
% ylabel('E (cm^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
% scrsz = get(0,'ScreenSize');
% set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
%     'Color',[1 1 1],'PaperPositionMode','auto')
% set(gca,'FontSize',13,'FontWeight','bold')
% title(['\epsilon = ',num2str(epsilon),' cm^2.s^{-3}'],'FontSize',11)
% 
% disp(['#--> Here comes some quantities of ',Name,' :'])
% disp(['epsilon = ',num2str(epsilon),' cm^2.s^{-3}'])
% disp(['L_beta = ',num2str(L_beta)])



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                          MODAL SPECTRA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Emn_mt = mean(reshape(Emn_t,Nt,Nrk,nTime),3);
% % % % 
% % % % % % -------------------------------- plots spectra
% % % % fig=figure; 
% % % % hax=axes; 
% % % % % loglog(n,(ERn(:,1)),'k','LineWidth',2)
% % % % loglog(n,(Emn_mt(1,:)),'-r','LineWidth',2) % I apply a factor 2, instead of doing the sum over m and -m
% % % % %
% % % % hold on
% % % % %
% % % % loglog(n,(2.*Emn_mt(2,:)),'-o','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
% % % % loglog(n,(2.*Emn_mt(3,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
% % % % loglog(n,(2.*Emn_mt(4,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
% % % % loglog(n,(2.*Emn_mt(5,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
% % % % loglog(n,(2.*Emn_mt(6,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
% % % % loglog(n,(2.*Emn_mt(7,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
% % % % loglog(n,(2.*Emn_mt(8,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
% % % % loglog(n,(2.*Emn_mt(9,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
% % % % loglog(n,(2.*Emn_mt(10,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
% % % % loglog(n,(2.*Emn_mt(11,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
% % % % legend('m=0','m = 1','m = 2','m = 3','m = 4','m = 5','m = 6','m = 7','m = 8','m = 9','m = 10')
% % % % % loglog(n,EZ_Theo,'--r')
% % % % % loglog(n,ER_Theo,'--k')
% % % % ylim([10.^(-7.) 0.1])
% % % % xlim([n(1) n(end)])
% % % % line([n_mag n_mag],get(hax,'YLim'),'Color',[0 0 0])
% % % % line([2*n_mag 2*n_mag],get(hax,'YLim'),'Color',[0 0 0])
% % % % line([n_jets n_jets],get(hax,'YLim'),'Color',[0 0 0])
% % % % box on
% % % % scrsz = get(0,'ScreenSize');
% % % % set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
% % % %     'Color',[1 1 1],'PaperPositionMode','auto')
% % % % set(gca,'FontSize',13,'FontWeight','bold')
% % % % xlabel('n','FontSize',13,'FontWeight','bold','Color','k')
% % % % ylabel('E','FontSize',13,'FontWeight','bold','Color','k')