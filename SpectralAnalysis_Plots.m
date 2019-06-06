clear all
close all
roots = '/media/simon/simon/ESP_29/'; % Root path..
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                             LOAD INFOS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run([roots,'InfosFile.m'])
load([roots,'/SpectralAnalysis_infos.mat'])
load('besselzeros2_C.mat'); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                              LOAD DATA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[EZn,a] = loadmtx([roots,'/Ezn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nFrames),'_window_',num2str(windowing)]);
[ERn,a] = loadmtx([roots,'/ERn_nbm',num2str(M),'_nbn',num2str(Nrk),'_Fr_',num2str(nFrames),'_window_',num2str(windowing)]);

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
% -------------------------------- plots spectra
figure
loglog(n,mean(EZn,2),'r')
hold on
loglog(n,mean(ERn,2),'k','LineWidth',2)
loglog(n,EZ_Theo,'--r')
loglog(n,ER_Theo,'--k')
ylim([10.^(-7.) 0.1])
xlim([n(1) n(end)])
box on
xlabel('n','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E (cm^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
set(gca,'FontSize',13,'FontWeight','bold')
title(['\epsilon = ',num2str(epsilon),' cm^2.s^{-3}'],'FontSize',11)