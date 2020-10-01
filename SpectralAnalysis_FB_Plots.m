load([roots,Name,'/SpectralAnalysis_infos.mat'])
% Itime=501;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              LOAD DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[EZn,a] = loadmtx([roots,Name,'/Ezn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nTime),'_Itime_',num2str(Itime)]);
[ERn,a] = loadmtx([roots,Name,'/ERn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nTime),'_Itime_',num2str(Itime)]);
if(0==0)
[Emn_t,a] = loadmtx([roots,Name,'/Emn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nTime),'_Itime_',num2str(Itime)]);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                SPECTRA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = [1:Nrk]; % non-dimensional
% cm2tom2 = 10.^(-4);
cm2tom2 = 1;
% -------------------------------- Theoretical spectra
EZ_Theo = (Cz*beta.^2.).*(1./R).*(J_root(1,1:Nrk)./(R)).^(-5.);
ER_Theo = Ck.*(epsilon.^(2./3.))*(1./R).*(J_root(2,1:Nrk)./(R)).^(-5./3.);
% ---------------------------------scales
L_beta=(epsilon./beta.^3.).^(1./5.);
Lo_beta = (Ck/Cz).^(3./10.).*(J_root(1,1:Nrk)./J_root(2,1:Nrk)).^(1./2.);
Lhat_beta = (J_root(1,1:Nrk)./n).*Lo_beta.*L_beta;
L_R = sqrt(2.*Urms./beta); 
n_r=R/L_R;
L_mag = 1.2 % magnet size in cm
L_jets = 7 % jets typical size in cm
n_mag=R/L_mag;
n_jets=R/L_jets;
% -------------------------------- plots spectra
fig=figure; 
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')
axes('FontSize',18,'Linewidth',2,'FontName','times',...
     'TickLength',[0.01; 0.03],'Xscale','Log',...
     'Position',[0.13 0.17 0.82 0.75])   
%
loglog(n,mean(EZn,2).*cm2tom2,'r')
hold on
loglog(n,mean(ERn,2).*cm2tom2,'k','LineWidth',2)
% loglog(n,mean(EZn+ERn,2).*cm2tom2,'k','LineWidth',2)
loglog(n,EZ_Theo.*cm2tom2,'--r')
loglog(n,ER_Theo.*cm2tom2,'--k')
ylim([10.^(-5.) 10.^3].*cm2tom2)
xlim([n(1) n(end)])
line([R/Lhat_beta(nhat) R/Lhat_beta(nhat)],get(gca,'YLim'),'Color',[0 0 0])
% line([n_r n_r],get(gca,'YLim'),'Color',[1 0 0])
% line([n_mag n_mag],get(gca,'YLim'),'Color',[0 0 0])
% line([2*n_mag 2*n_mag],get(gca,'YLim'),'Color',[0 0 0])
box on
xlabel ('$n$','FontSize',18,'FontName','times','Interpreter','Latex')
ylabel ('$E $(m$^{2}$s$^{-2}$)','FontSize',18,'FontName','times','Interpreter','Latex')

title(['$\varepsilon =$ ',num2str(epsilon),' (m$^{2}$s$^{-3}$)'],'FontSize',18,'FontName','times','Interpreter','Latex')

disp(['#--> Here comes some quantities of ',Name,' :'])
disp(['epsilon = ',num2str(epsilon),' cm^2.s^{-3}'])
disp(['L_beta = ',num2str(L_beta)])
disp(['Lhat_beta = ',num2str(Lhat_beta(nhat))])
disp(['Con n = ',num2str(nhat)])
disp(['L_R = ',num2str(L_R)])
disp(['TKE = ',num2str(sum(mean(EZn(2:end,:)+ERn(2:end,:),2).*cm2tom2)),' m^2.s^{-2}'])


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                          MODAL SPECTRA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(0==0)
Emn_mt = mean(reshape(Emn_t,Nt,Nrk,nTime),3);
% Emn = reshape(Emn_t,Nt,Nrk,nTime);
% Emn_mt = mean(Emn(:,:,800:1501),3);

% % -------------------------------- plots spectra
fig=figure; 
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')
axes('FontSize',18,'Linewidth',2,'FontName','times',...
     'TickLength',[0.01; 0.03],'Xscale','Log',...
     'Position',[0.13 0.17 0.82 0.75])   
%
loglog(n,(Emn_mt(1,:)),'-r','LineWidth',2) % I apply a factor 2, instead of doing the sum over m and -m
%
hold on
%
loglog(n,(2.*Emn_mt(2,:)),'-o','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(3,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(4,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(5,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(6,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(7,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(8,:)),'->','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(9,:)),'-*','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(10,:)),'-<','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn_mt(11,:)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
legend('m=0','m = 1','m = 2','m = 3','m = 4','m = 5','m = 6','m = 7','m = 8','m = 9','m = 10')
% loglog(n,EZ_Theo,'--r')
% loglog(n,ER_Theo,'--k')
% ylim([10.^(-7.) 0.1])
xlim([n(1) n(end)])
line([n_mag n_mag],get(gca,'YLim'),'Color',[0 0 0])
line([2*n_mag 2*n_mag],get(gca,'YLim'),'Color',[0 0 0])
line([n_jets n_jets],get(gca,'YLim'),'Color',[0 0 0])
box on
%
xlabel ('$n$','FontSize',18,'FontName','times','Interpreter','Latex')
ylabel ('$E $(m$^{2}$s$^{-2}$)','FontSize',18,'FontName','times','Interpreter','Latex')
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                MODAL SPECTRA temporale:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(1==0)
Emn = reshape(Emn_t,Nt,Nrk,nTime);
tt=490;
% % -------------------------------- plots spectra
fig=figure; 
hax=axes; 
% loglog(n,(ERn(:,1)),'k','LineWidth',2)
loglog(n,(Emn(1,:,tt))+(2.*Emn(2,:,tt)),'-r','LineWidth',2) % I apply a factor 2, instead of doing the sum over m and -m
%
hold on
%
% loglog(n,(2.*Emn(2,:,tt)),'-k','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn(3,:,tt)),'-m','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn(4,:,tt)),'-k','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn(5,:,tt)),'-k','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn(6,:,tt)),'-k','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn(7,:,tt)),'-og','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn(8,:,tt)),'->c','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
loglog(n,(2.*Emn(9,:,tt)),'-*b','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
% loglog(n,(2.*Emn(10,:,tt)),'-<','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
% loglog(n,(2.*Emn(11,:,tt)),'-','LineWidth',1) % I apply a factor 2, instead of doing the sum over m and -m
% legend('m=0','m = 1','m = 2','m = 3','m = 4','m = 5','m = 6','m = 7','m = 8')
legend('m=0 + m = 1','m = 2','m = 3','m = 4','m = 5','m = 6','m = 7','m = 8')
% loglog(n,EZ_Theo,'--r')
% loglog(n,ER_Theo,'--k')
% ylim([10.^(-7.) 0.1])
xlim([n(1) n(end)])
% line([n_mag n_mag],get(hax,'YLim'),'Color',[0 0 0])
% line([2*n_mag 2*n_mag],get(hax,'YLim'),'Color',[0 0 0])
% line([n_jets n_jets],get(hax,'YLim'),'Color',[0 0 0])
box on
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')
set(gca,'FontSize',13,'FontWeight','bold')
xlabel('n','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E (m^2.s^{-2})','FontSize',13,'FontWeight','bold','Color','k')
end