load([roots,Name,'/EFluxes_Spectral_infos.mat'])
addpath('/home/simon/Bureau/Esperimento-DICEA/JUMP/Altre-funzioni/cbrewer/')  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              LOAD DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Energy_flux_mt_yxl,a] = loadmtx([roots,Name,'/Energy_flux_maps_Fr_',num2str(nTime),'_Itime_',num2str(Itime),'_Idymin_',num2str(idymin)]);
%
cmax = 5.e-06;
cmax = 5.e-04;
%--------------------------------------------------
figure
% contourf(X(Ay,Ax),Y(Ay,Ax),FLUSSO_Energia_tc(:,:,4), 250, 'LineStyle','none');
ids = 2;
imagesc(FLUSSO_Energia_tc(:,:,ids));
title(['Energy flux at scale l = ',num2str(ids)],'FontSize',18,'FontWeight','bold','Interpreter','Latex')
h = colorbar;
set(get(h,'label'),'string','$\Pi_{\varepsilon}$(m$^2$s$^{-3}$)','FontSize',18,'FontWeight','bold','Interpreter','Latex');
caxis([-cmax cmax]);
%--------------------------------------------------
figure
% contourf(X(Ay,Ax),Y(Ay,Ax),FLUSSO_Energia_tc(:,:,4), 250, 'LineStyle','none');
ids = 3;
imagesc(FLUSSO_Energia_tc(:,:,ids));
title(['Energy flux at scale l = ',num2str(ids)],'FontSize',18,'FontWeight','bold','Interpreter','Latex')
h = colorbar;
set(get(h,'label'),'string','$\Pi_{\varepsilon}$(m$^2$s$^{-3}$)','FontSize',18,'FontWeight','bold','Interpreter','Latex');
caxis([-cmax cmax]);
%--------------------------------------------------
figure
% contourf(X(Ay,Ax),Y(Ay,Ax),FLUSSO_Energia_tc(:,:,4), 250, 'LineStyle','none');
ids = 6;
imagesc(FLUSSO_Energia_tc(:,:,ids));
title(['Energy flux at scale l = ',num2str(ids)],'FontSize',18,'FontWeight','bold','Interpreter','Latex')
h = colorbar;
set(get(h,'label'),'string','$\Pi_{\varepsilon}$(m$^2$s$^{-3}$)','FontSize',18,'FontWeight','bold','Interpreter','Latex');
caxis([-cmax cmax]);
%--------------------------------------------------
figure
% contourf(X(Ay,Ax),Y(Ay,Ax),FLUSSO_Energia_tc(:,:,4), 250, 'LineStyle','none');
ids = 20;
imagesc(FLUSSO_Energia_tc(:,:,ids));
title(['Energy flux at scale l = ',num2str(ids)],'FontSize',18,'FontWeight','bold','Interpreter','Latex')
h = colorbar;
set(get(h,'label'),'string','$\Pi_{\varepsilon}$(m$^2$s$^{-3}$)','FontSize',18,'FontWeight','bold','Interpreter','Latex');
caxis([-cmax cmax]);

% 
% figure
% contourf(X(Ay,Ax),Y(Ay,Ax),FLUSSO_Energia_tc(:,:,18), 250, 'LineStyle','none');
% title('Energy flux - large scale  in the averaged surface')
% colorbar
% caxis([-0.01 0.01]);