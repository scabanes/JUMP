load([roots,Name,'/EFluxes_Spectral_infos.mat'])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              LOAD DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Fluxes_mt,a] = loadmtx([roots,Name,'/Fluxes_mt']);
[Energy_flux_t,a] = loadmtx([roots,Name,'/Energy_flux_Fr_',num2str(nTime)]);
[Enstrophy_flux_t,a] = loadmtx([roots,Name,'/Enstrophy_flux_Fr_',num2str(nTime)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           EFluces.plot
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
% semilogx(2.*pi./(vl.*dx),FLUSSO_Energia_tl,'c')
semilogx(2.*pi./(vl.*dx),Fluxes_mt(:,1),'-b','linewidth',3)
% semilogx((2.*pi.*[Nkx-1:-0.5:0])/Lx,FLUSSO_Energia_l,'LineWidth',3,'Color','k')
hold on
semilogx(2.*pi./(vl.*dx),Fluxes_mt(:,2)./(10.^(2)),'-k','linewidth',3)
% semilogx((2.*pi.*[Nkx-1:-0.5:0])/Lx,FLUSSO_Enstrofia_l./(10.^(2)),'LineWidth',3,'Color','r')
yline(0,':');
grid on
box on
xlim([2.*pi./(vl(end).*dx) 2.*pi./(vl(1).*dx)])
%----- scales
% forcing, magnets
% % line([1.2 1.2],get(gca,'YLim'),'Color','k')
% % line([2*1.2 2*1.2],get(gca,'YLim'),'Color','k')
% % % jets/Vortices separation
% % line([6 6],get(gca,'YLim'),'Color','k')
ylabel ('$\Pi$(m$^2$s$^{-3}$), $Z.10^2$(s$^{-3}$)','FontSize',18,'FontName','times','Interpreter','Latex')
%xlabel ('$k$(cm$^{-1}$)','FontSize',18,'FontName','times','Interpreter','Latex')
xlabel ('$k$(m$^{-1}$)','FontSize',18,'FontName','times','Interpreter','Latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Time evolution
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
semilogx(2.*pi./(vl.*dx),Energy_flux_t,'c')
hold on
semilogx(2.*pi./(vl.*dx),Fluxes_mt(:,1),'LineWidth',3,'Color','b')
yline(0,':');
grid on
box on
xlim([2.*pi./(vl(end).*dx) 2.*pi./(vl(3).*dx)])
%----- scales
% forcing, magnets
% % line([1.2 1.2],get(gca,'YLim'),'Color','k')
% % line([2*1.2 2*1.2],get(gca,'YLim'),'Color','k')
% % % jets/Vortices separation
% % line([6 6],get(gca,'YLim'),'Color','k')
ylabel ('$\Pi$(m$^2$s$^{-3}$), $Z$(s$^{-3}$)','FontSize',18,'FontName','times','Interpreter','Latex')
%xlabel ('$k$(cm$^{-1}$)','FontSize',18,'FontName','times','Interpreter','Latex')
xlabel ('$k$(m$^{-1}$)','FontSize',18,'FontName','times','Interpreter','Latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

