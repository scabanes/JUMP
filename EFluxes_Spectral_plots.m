load([roots,Name,'/EFluxes_Spectral_infos.mat'])

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              LOAD DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Fluxes_mt,a] = loadmtx([roots,Name,'/Fluxes_mt_Fr_',num2str(nTime),'_Itime_',num2str(Itime),'_Idymin_',num2str(idymin)]);
[Energy_flux_t,a] = loadmtx([roots,Name,'/Energy_flux_Fr_',num2str(nTime),'_Itime_',num2str(Itime),'_Idymin_',num2str(idymin)]);
[Enstrophy_flux_t,a] = loadmtx([roots,Name,'/Enstrophy_flux_Fr_',num2str(nTime),'_Itime_',num2str(Itime),'_Idymin_',num2str(idymin)]);

factEnstro = 1;% 10.^(3);%10.^(-12);
cmTom = 1;%10.^(-4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           EFluces.plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% scrsz = get(0,'ScreenSize');
% set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
%     'Color',[1 1 1],'PaperPositionMode','auto')
% axes('FontSize',18,'Linewidth',2,'FontName','times',...
%      'TickLength',[0.01; 0.03],'Xscale','Log',...
%      'Position',[0.13 0.17 0.82 0.75])   
hold on
subplot ('Position',[0.13 0.17 0.82 0.75]) 
% semilogx(2.*pi./(vl.*dx),FLUSSO_Energia_tl,'c')
semilogx(Nx./vl,Fluxes_mt(:,1).*cmTom,'-b','linewidth',3)
% semilogx(Nx./vl,Energy_flux_t.*cmTom,'-b','linewidth',3)
% semilogx((2.*pi.*[Nkx-1:-0.5:0])/Lx,FLUSSO_Energia_l,'LineWidth',3,'Color','k')
hold on
semilogx(Nx./vl,Fluxes_mt(:,2).*factEnstro,'-k','linewidth',3)
% semilogx(Nx./vl,Enstrophy_flux_t.*factEnstro,'-k','linewidth',3)
% semilogx((2.*pi.*[Nkx-1:-0.5:0])/Lx,FLUSSO_Enstrofia_l./(10.^(2)),'LineWidth',3,'Color','r')
yline(0,':','linewidth',1);
% grid on
box on
% xlim([Nx/vl(1) Nx/vl(end)])
%----- scales
xlabel('n','FontSize',13,'FontWeight','bold','Color','k')
ylabel (['$\Pi_{\epsilon}$(m$^2$s$^{-3}$), $Z \times ',num2str(factEnstro),'$(s$^{-3}$)'],'FontSize',18,'FontWeight','bold','Interpreter','Latex')
scrsz = get(0,'ScreenSize');
% set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
%     'Color',[1 1 1],'PaperPositionMode','auto')
set(gca,'FontSize',13,'FontWeight','bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Time evolution Energy flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')
axes('FontSize',18,'Linewidth',2,'FontName','times',...
     'TickLength',[0.01; 0.03],'Xscale','Log',...
     'Position',[0.13 0.17 0.82 0.75])   
hold on
subplot ('Position',[0.13 0.17 0.82 0.75]) 
semilogx(1./(vl.*dx),Energy_flux_t.*cmTom,'c')
hold on
% semilogx(Nx./vl,Fluxes_mt(:,1).*cmTom,'LineWidth',3,'Color','b')
semilogx(1./(vl.*dx),Fluxes_mt(:,1).*cmTom,'LineWidth',3,'Color','b')
% semilogx(vl,Fluxes_mt(:,1).*cmTom,'LineWidth',3,'Color','b')
% hold on
yline(0,':');
% grid on
box on
% xlim([(2.*pi.*Nx*dx)./(vl(end).*dx) (2.*pi.*Nx*dx)./(vl(1).*dx)])
xlim([1./(vl(end).*dx) 1./(vl(1).*dx)])
%----- scales
% forcing, magnets
% % line([1.2 1.2],get(gca,'YLim'),'Color','k')
% % line([2*1.2 2*1.2],get(gca,'YLim'),'Color','k')
% % % jets/Vortices separation
% % line([6 6],get(gca,'YLim'),'Color','k')
line([1./(2.*dx) 1./(2.*dx)],get(gca,'YLim'),'Color','k')
line([1./(3.*dx) 1./(3.*dx)],get(gca,'YLim'),'Color','k')
line([1./(6.*dx) 1./(6.*dx)],get(gca,'YLim'),'Color','k')
ylabel ('$\Pi$(m$^2$s$^{-3}$)','FontSize',18,'FontName','times','Interpreter','Latex')
%xlabel ('$k$(cm$^{-1}$)','FontSize',18,'FontName','times','Interpreter','Latex')
% xlabel ('$n$','FontSize',18,'FontName','times','Interpreter','Latex')
xlabel ('$k$ [m$^{-1}$]','FontSize',18,'FontName','times','Interpreter','Latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Time evolution Enstrophy Flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
scrsz = get(0,'ScreenSize');
set(gcf,'Position',[0 scrsz(4)/3 scrsz(3)/2.7 scrsz(4)/2.5],...
    'Color',[1 1 1],'PaperPositionMode','auto')
axes('FontSize',18,'Linewidth',2,'FontName','times',...
     'TickLength',[0.01; 0.03],'Xscale','Log',...
     'Position',[0.13 0.17 0.82 0.75])   
hold on
subplot ('Position',[0.13 0.17 0.82 0.75]) 
semilogx(1./(vl.*dx),Enstrophy_flux_t,'c')
hold on
% semilogx(Nx./vl,Fluxes_mt(:,2),'LineWidth',3,'Color','k') % addimentional wavenumber
semilogx(1./(vl.*dx),Fluxes_mt(:,2),'LineWidth',3,'Color','k') % enstrophy = f(m-1)
% hold on
yline(0,':');
% grid on
box on
% xlim([(2.*pi.*Nx*dx)./(vl(end).*dx) (2.*pi.*Nx*dx)./(vl(1).*dx)])
xlim([1./(vl(end).*dx) 1./(vl(1).*dx)])
%----- scales
% forcing, magnets
% % line([1.2 1.2],get(gca,'YLim'),'Color','k')
% % line([2*1.2 2*1.2],get(gca,'YLim'),'Color','k')
% % % jets/Vortices separation
% % line([6 6],get(gca,'YLim'),'Color','k')
ylabel ('$Z $(s$^{-3}$)','FontSize',18,'FontName','times','Interpreter','Latex')
%xlabel ('$k$(cm$^{-1}$)','FontSize',18,'FontName','times','Interpreter','Latex')
% % % xlabel ('$n$','FontSize',18,'FontName','times','Interpreter','Latex')
xlabel ('$k$ [m$^{-1}$]','FontSize',18,'FontName','times','Interpreter','Latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
