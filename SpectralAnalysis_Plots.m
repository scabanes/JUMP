load([roots,Name,'/SpectralAnalysis_infos.mat'])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              LOAD DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[EZn,a] = loadmtx([roots,Name,'/Ezn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
[ERn,a] = loadmtx([roots,Name,'/ERn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nTime),'_window_',num2str(windowing)]);
[Emn_t,a] = loadmtx([roots,Name,'/Emn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nTime),'_window_',num2str(windowing)]);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                SPECTRA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = [1:Nrk]; % non-dimensional
% -------------------------------- Theoretical spectra
EZ_Theo = (Cz*beta.^2.).*(1./R).*(J_root(1,1:Nrk)./R).^(-5.);
ER_Theo = Ck.*(epsilon.^(2./3.))*(1./R).*(J_root(2,1:Nrk)./R).^(-5./3.);
% ---------------------------------scales
L_beta=(epsilon./beta.^3.).^(1./5.);
Lo_beta = (Ck/Cz).^(3./10.).*(J_root(1,1:Nrk)./J_root(2,1:Nrk)).^(1./2.);
Lhat_beta = (J_root(1,1:Nrk)./n).*Lo_beta.*L_beta;
L_R = sqrt(2.*Urms./beta); 
n_r=R/L_R;
L_mag = 1.2 % magnet size in cm
n_mag=R/L_mag;
% -------------------------------- plots spectra
fig=figure; 
hax=axes; 
loglog(n,mean(EZn,2),'r')
hold on
loglog(n,mean(ERn,2),'k','LineWidth',2)
loglog(n,EZ_Theo,'--r')
loglog(n,ER_Theo,'--k')
ylim([10.^(-7.) 0.3])
xlim([n(1) n(end)])
%line([R/Lhat_beta(nint) R/Lhat_beta(nint)],get(hax,'YLim'),'Color',[0 0 0])
% line([n_r n_r],get(hax,'YLim'),'Color',[0 0 0])
line([n_mag n_mag],get(hax,'YLim'),'Color',[0 0 0])
line([2*n_mag 2*n_mag],get(hax,'YLim'),'Color',[0 0 0])
box on
xlabel('n','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E (cm^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')
set(gca,'FontSize',13,'FontWeight','bold')
title(['\epsilon = ',num2str(epsilon),' cm^2.s^{-3}'],'FontSize',11)

disp(['#--> Here comes some quantities of ',Name,' :'])
disp(['epsilon = ',num2str(epsilon),' cm^2.s^{-3}'])
disp(['L_beta = ',num2str(L_beta)])
disp(['Lhat_beta = ',num2str(Lhat_beta(nhat))])
disp(['Con n = ',num2str(nhat)])
disp(['L_R = ',num2str(L_R)])


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                          MODAL SPECTRA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Emn_mt = mean(reshape(Emn_t,Nt,Nrk,nTime),3);

% % -------------------------------- plots spectra
fig=figure; 
hax=axes; 
loglog(n,(ERn(:,1)),'k','LineWidth',2)
hold on
loglog(n,(Emn_mt(1,:)),'-r','LineWidth',2) % I apply a factor 2, instead of doing the sum over m and -m
%
loglog(n,(2.*Emn_mt(2,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(3,:)),':','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(4,:)),':','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(5,:)),':','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(6,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(7,:)),':','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(8,:)),':','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(9,:)),':','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(10,:)),':','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(11,:)),':','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
legend('ER','m=0','m = 1','m = 2','m = 3','m = 4','m = 5','m = 6','m = 7','m = 8','m = 9','m = 10')
% loglog(n,EZ_Theo,'--r')
% loglog(n,ER_Theo,'--k')
ylim([10.^(-7.) 0.3])
xlim([n(1) n(end)])
  scrsz = get(0,'ScreenSize');

   set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')
set(gca,'FontSize',13,'FontWeight','bold')
