clear all
close all

%**********************************************************
%           PARAMETRI DA SCEGLIERE
%**********************************************************
expt = {'5';'6';'7';'8';'9';'3'}
roots = '../Dati-Torino/';
output = 'outputs'
windowing=1;
load('D:\simon\script\besselzeros2_C.mat'); 
EXPT_infos
LXAXIS=[10.^(-7.) 0.04] % X axis limits
%**********************************************************
L_beta_exp = zeros(length(expt),1);
Lhat_beta_exp = zeros(length(expt),1);
LM_exp = zeros(length(expt),1);
n_intersec = zeros(length(expt),1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                  EXP 3:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------
%------------------------------ windowing
cexp = 6;
load([roots,names{cexp,2},'/',output,'/infos.mat'])
nFrames = ((endi{cexp,2}-starti{cexp,2})/step + 1);
load([roots,names{cexp,2},'/',output,'/E_infos.mat'])
R=R_trunc;% in centimetri
%%%%%%%%%%
name_mtx1 = [roots,names{cexp,2},'/',output,'/Ezn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(nFrames),'_window',num2str(windowing)];
[Ezn,b]=loadmtx(name_mtx1);
name_mtx2 = [roots,names{cexp,2},'/',output,'/ERn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(nFrames),'_window',num2str(windowing)];
[ERn,b]=loadmtx(name_mtx2);

Cz = 0.5;
Ck = 6;
beta_exp = beta{cexp,2}/100.; % cm-1.s-1
epsilon = 1.5*10.^(-5); % cm-2.s-3
EZ_Theo = (Cz*beta_exp.^2.).*(1./R).*(J_root(1,1:Nrk)./R).^(-5.);
ER_Theo = Ck.*(epsilon.^(2./3.))*(1./R).*(J_root(2,1:Nrk)./R).^(-5./3.);
U=0.0044*100;%cm
L_beta=(epsilon./beta_exp.^3.).^(1./5.);
Lo_beta = (Ck/Cz).^(3./10.).*(J_root(1,1:Nrk)./J_root(2,1:Nrk)).^(1./2.).*L_beta;
Lhat_beta = (J_root(1,1:Nrk)./n).*Lo_beta;
L_R = sqrt(2.*U./beta_exp); 
n_r=R/L_R;
%
n_intersec(cexp)=5;
L_beta_exp(cexp) = L_beta;
Lhat_beta_exp(cexp) = Lhat_beta(n_intersec(cexp));
LM_exp(cexp) = 0.072*100;% in centimetri
% -------------------------------- plots spectra
figure;hold on
plot(n,mean(Ezn,2),'r','linewidth',2.)
hold on
plot(n,mean(ERn,2),'k','linewidth',2.)
plot(n,EZ_Theo,'--r')
plot(n,ER_Theo,'--k')
% line([10 10], get(gca, 'ylim'));
set(gca,'YScale', 'log','XScale', 'log')
ylim(LXAXIS)
xlim([n(1) n(end)])
box on
xlabel('n','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E (cm^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
set(gca,'FontSize',13,'FontWeight','bold')
title([names{cexp,2},': \epsilon = ',num2str(epsilon),' cm^2.s^{-3}, windowing = ',num2str(windowing)],'FontSize',11)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                  EXP 9:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------
%------------------------------ windowing
cexp = 5;
load([roots,names{cexp,2},'/',output,'/infos.mat'])
nFrames = ((endi{cexp,2}-starti{cexp,2})/step + 1);
load([roots,names{cexp,2},'/',output,'/E_infos.mat'])
R=R_trunc;
%%%%%%%%%%
name_mtx1 = [roots,names{cexp,2},'/',output,'/Ezn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(nFrames),'_window',num2str(windowing)];
[Ezn,b]=loadmtx(name_mtx1);
name_mtx2 = [roots,names{cexp,2},'/',output,'/ERn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(nFrames),'_window',num2str(windowing)];
[ERn,b]=loadmtx(name_mtx2);

Cz = 0.5;
Ck = 6;
beta_exp = beta{cexp,2}/100.; % cm-1.s-1
epsilon = 2.5*10.^(-5); % cm-2.s-3
EZ_Theo = (Cz*beta_exp.^2.).*(1./R).*(J_root(1,1:Nrk)./R).^(-5.);
ER_Theo = Ck.*(epsilon.^(2./3.))*(1./R).*(J_root(2,1:Nrk)./R).^(-5./3.);
U=0.0044*100;%cm
L_beta=(epsilon./beta_exp.^3.).^(1./5.);
Lo_beta = (Ck/Cz).^(3./10.).*(J_root(1,1:Nrk)./J_root(2,1:Nrk)).^(1./2.).*L_beta;
Lhat_beta = (J_root(1,1:Nrk)./n).*Lo_beta;
L_R = sqrt(2.*U./beta_exp); 
n_r=R/L_R;
%
n_intersec(cexp)=5;
L_beta_exp(cexp) = L_beta;
Lhat_beta_exp(cexp) = Lhat_beta(n_intersec(cexp));
LM_exp(cexp) = 0.0971*100;% in centimetri
% -------------------------------- plots spectra
figure;hold on
plot(n,mean(Ezn,2),'r','linewidth',2.)
hold on
plot(n,mean(ERn,2),'k','linewidth',2.)
plot(n,EZ_Theo,'--r')
plot(n,ER_Theo,'--k')
% line([10 10], get(gca, 'ylim'));
set(gca,'YScale', 'log','XScale', 'log')
ylim(LXAXIS)
xlim([n(1) n(end)])
box on
xlabel('n','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E (cm^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
set(gca,'FontSize',13,'FontWeight','bold')
title([names{cexp,2},': \epsilon = ',num2str(epsilon),' cm^2.s^{-3}, windowing = ',num2str(windowing)],'FontSize',11)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                  EXP 7:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------
%------------------------------ windowing
cexp = 3;
load([roots,names{cexp,2},'/',output,'/infos.mat'])
nFrames = ((endi{cexp,2}-starti{cexp,2})/step + 1);
load([roots,names{cexp,2},'/',output,'/E_infos.mat'])
R=R_trunc;
%%%%%%%%%%
name_mtx1 = [roots,names{cexp,2},'/',output,'/Ezn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(nFrames),'_window',num2str(windowing)];
[Ezn,b]=loadmtx(name_mtx1);
name_mtx2 = [roots,names{cexp,2},'/',output,'/ERn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(nFrames),'_window',num2str(windowing)];
[ERn,b]=loadmtx(name_mtx2);

Cz = 0.5;
Ck = 6;
beta_exp = beta{cexp,2}/100.; % cm-1.s-1
epsilon = 1.65*10.^(-5); % cm-2.s-3
EZ_Theo = (Cz*beta_exp.^2.).*(1./R).*(J_root(1,1:Nrk)./R).^(-5.);
ER_Theo = Ck.*(epsilon.^(2./3.))*(1./R).*(J_root(2,1:Nrk)./R).^(-5./3.);
U=0.0044*100;%cm
L_beta=(epsilon./beta_exp.^3.).^(1./5.);
Lo_beta = (Ck/Cz).^(3./10.).*(J_root(1,1:Nrk)./J_root(2,1:Nrk)).^(1./2.).*L_beta;
Lhat_beta = (J_root(1,1:Nrk)./n).*Lo_beta;
L_R = sqrt(2.*U./beta_exp); 
n_r=R/L_R;
%
n_intersec(cexp)=4;
L_beta_exp(cexp) = L_beta;
Lhat_beta_exp(cexp) = Lhat_beta(n_intersec(cexp));
LM_exp(cexp) = 0.101*100;% in centimetri
% -------------------------------- plots spectra
figure;hold on
plot(n,mean(Ezn,2),'r','linewidth',2.)
hold on
plot(n,mean(ERn,2),'k','linewidth',2.)
plot(n,EZ_Theo,'--r')
plot(n,ER_Theo,'--k')
% line([10 10], get(gca, 'ylim'));
set(gca,'YScale', 'log','XScale', 'log')
ylim(LXAXIS)
xlim([n(1) n(end)])
box on
xlabel('n','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E (cm^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
set(gca,'FontSize',13,'FontWeight','bold')
title([names{cexp,2},': \epsilon = ',num2str(epsilon),' cm^2.s^{-3}, windowing = ',num2str(windowing)],'FontSize',11)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                  EXP 5:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------
%------------------------------ windowing
cexp = 1;
load([roots,names{cexp,2},'/',output,'/infos.mat'])
nFrames = ((endi{cexp,2}-starti{cexp,2})/step + 1);
load([roots,names{cexp,2},'/',output,'/E_infos.mat'])
R=R_trunc;
%%%%%%%%%%
name_mtx1 = [roots,names{cexp,2},'/',output,'/Ezn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(nFrames),'_window',num2str(windowing)];
[Ezn,b]=loadmtx(name_mtx1);
name_mtx2 = [roots,names{cexp,2},'/',output,'/ERn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(nFrames),'_window',num2str(windowing)];
[ERn,b]=loadmtx(name_mtx2);

Cz = 0.5;
Ck = 6;
beta_exp = beta{cexp,2}/100.; % cm-1.s-1
epsilon = 0.8*10.^(-5); % cm-2.s-3
EZ_Theo = (Cz*beta_exp.^2.).*(1./R).*(J_root(1,1:Nrk)./R).^(-5.);
ER_Theo = Ck.*(epsilon.^(2./3.))*(1./R).*(J_root(2,1:Nrk)./R).^(-5./3.);
U=0.0044*100;%cm
L_beta=(epsilon./beta_exp.^3.).^(1./5.);
Lo_beta = (Ck/Cz).^(3./10.).*(J_root(1,1:Nrk)./J_root(2,1:Nrk)).^(1./2.).*L_beta;
Lhat_beta = (J_root(1,1:Nrk)./n).*Lo_beta;
L_R = sqrt(2.*U./beta_exp); 
n_r=R/L_R;
%
n_intersec(cexp)=3;
L_beta_exp(cexp) = L_beta;
Lhat_beta_exp(cexp) = Lhat_beta(n_intersec(cexp));
LM_exp(cexp) = 0.17*100;% in centimetri
% -------------------------------- plots spectra
figure;hold on
plot(n,mean(Ezn,2),'r','linewidth',2.)
hold on
plot(n,mean(ERn,2),'k','linewidth',2.)
plot(n,EZ_Theo,'--r')
plot(n,ER_Theo,'--k')
% line([10 10], get(gca, 'ylim'));
set(gca,'YScale', 'log','XScale', 'log')
ylim(LXAXIS)
xlim([n(1) n(end)])
box on
xlabel('n','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E (cm^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
set(gca,'FontSize',13,'FontWeight','bold')
title([names{cexp,2},': \epsilon = ',num2str(epsilon),' cm^2.s^{-3}, windowing = ',num2str(windowing)],'FontSize',11)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                  EXP 8:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------
%------------------------------ windowing
cexp = 4;
load([roots,names{cexp,2},'/',output,'/infos.mat'])
nFrames = ((endi{cexp,2}-starti{cexp,2})/step + 1);
load([roots,names{cexp,2},'/',output,'/E_infos.mat'])
R=R_trunc;
%%%%%%%%%%
name_mtx1 = [roots,names{cexp,2},'/',output,'/Ezn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(2150),'_window',num2str(windowing)];
[Ezn,b]=loadmtx(name_mtx1);
name_mtx2 = [roots,names{cexp,2},'/',output,'/ERn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(2150),'_window',num2str(windowing)];
[ERn,b]=loadmtx(name_mtx2);

Cz = 0.5;
Ck = 6;
beta_exp = beta{cexp,2}/100.; % cm-1.s-1
epsilon = 2.*10.^(-5); % cm-2.s-3
EZ_Theo = (Cz*beta_exp.^2.).*(1./R).*(J_root(1,1:Nrk)./R).^(-5.);
ER_Theo = Ck.*(epsilon.^(2./3.))*(1./R).*(J_root(2,1:Nrk)./R).^(-5./3.);
U=0.0044*100;%cm
L_beta=(epsilon./beta_exp.^3.).^(1./5.);
Lo_beta = (Ck/Cz).^(3./10.).*(J_root(1,1:Nrk)./J_root(2,1:Nrk)).^(1./2.).*L_beta;
Lhat_beta = (J_root(1,1:Nrk)./n).*Lo_beta;
L_R = sqrt(2.*U./beta_exp); 
n_r=R/L_R;
%
n_intersec(cexp)=3;
L_beta_exp(cexp) = L_beta;
Lhat_beta_exp(cexp) = Lhat_beta(n_intersec(cexp));
LM_exp(cexp) = 0.161*100;% in centimetri
% -------------------------------- plots spectra
figure;hold on
plot(n,mean(Ezn,2),'r','linewidth',2.)
hold on
plot(n,mean(ERn,2),'k','linewidth',2.)
plot(n,EZ_Theo,'--r')
plot(n,ER_Theo,'--k')
% line([10 10], get(gca, 'ylim'));
set(gca,'YScale', 'log','XScale', 'log')
ylim(LXAXIS)
xlim([n(1) n(end)])
box on
xlabel('n','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E (cm^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
set(gca,'FontSize',13,'FontWeight','bold')
title([names{cexp,2},': \epsilon = ',num2str(epsilon),' cm^2.s^{-3}, windowing = ',num2str(windowing)],'FontSize',11)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                  EXP 6:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------
%------------------------------ windowing
cexp = 2;
load([roots,names{cexp,2},'/',output,'/infos.mat'])
nFrames = ((endi{cexp,2}-starti{cexp,2})/step + 1);
load([roots,names{cexp,2},'/',output,'/E_infos.mat'])
R=R_trunc;
%%%%%%%%%%
name_mtx1 = [roots,names{cexp,2},'/',output,'/Ezn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(nFrames),'_window',num2str(windowing)];
[Ezn,b]=loadmtx(name_mtx1);
name_mtx2 = [roots,names{cexp,2},'/',output,'/ERn_nbm',num2str(M),'_nbn',num2str(Nrk),'_',num2str(nFrames),'_window',num2str(windowing)];
[ERn,b]=loadmtx(name_mtx2);

Cz = 0.5;
Ck = 6;
beta_exp = beta{cexp,2}/100.; % cm-1.s-1
epsilon = 2.*10.^(-5); % cm-2.s-3
EZ_Theo = (Cz*beta_exp.^2.).*(1./R).*(J_root(1,1:Nrk)./R).^(-5.);
ER_Theo = Ck.*(epsilon.^(2./3.))*(1./R).*(J_root(2,1:Nrk)./R).^(-5./3.);
U=0.0044*100;%cm
L_beta=(epsilon./beta_exp.^3.).^(1./5.);
Lo_beta = (Ck/Cz).^(3./10.).*(J_root(1,1:Nrk)./J_root(2,1:Nrk)).^(1./2.).*L_beta;
Lhat_beta = (J_root(1,1:Nrk)./n).*Lo_beta;
L_R = sqrt(2.*U./beta_exp); 
n_r=R/L_R;
%
n_intersec(cexp)=2;
L_beta_exp(cexp) = L_beta;
Lhat_beta_exp(cexp) = Lhat_beta(n_intersec(cexp));
LM_exp(cexp) = 0.164*100;% in centimetri
% -------------------------------- plots spectra
figure;hold on
plot(n,mean(Ezn,2),'r','linewidth',2.)
hold on
plot(n,mean(ERn,2),'k','linewidth',2.)
plot(n,EZ_Theo,'--r')
plot(n,ER_Theo,'--k')
% line([10 10], get(gca, 'ylim'));
set(gca,'YScale', 'log','XScale', 'log')
ylim(LXAXIS)
xlim([n(1) n(end)])
box on
xlabel('n','FontSize',13,'FontWeight','bold','Color','k')
ylabel('E (cm^2.s^2)','FontSize',13,'FontWeight','bold','Color','k')
set(gca,'FontSize',13,'FontWeight','bold')
title([names{cexp,2},': \epsilon = ',num2str(epsilon),' cm^2.s^{-3}, windowing = ',num2str(windowing)],'FontSize',11)



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                        SYNTHESIS PLOTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA = LM_exp\ Lhat_beta_exp
figure;hold on
plot(Lhat_beta_exp, LM_exp, '>k','MarkerFaceColor','k')
plot([30:1:75],[30:1:75].*(1./AA),'k')
% legend('','y=ax')
xlabel('hat L_\beta','FontSize',13,'FontWeight','bold','Color','k')
ylabel('L_M','FontSize',13,'FontWeight','bold','Color','k')
set(gca,'FontSize',13,'FontWeight','bold')
%
BB = LM_exp\ L_beta_exp;
figure;hold on
plot(L_beta_exp,LM_exp, '>k','MarkerFaceColor','k')
plot([5:1:13],[5:1:13].*(1./BB),'k')
xlabel('L_\beta','FontSize',13,'FontWeight','bold','Color','k')
ylabel('L_M','FontSize',13,'FontWeight','bold','Color','k')
set(gca,'FontSize',13,'FontWeight','bold')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                  Print:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% SCALES
disp('############# SCAlES ################################')
disp('LM = ')
disp(num2str(LM_exp))
disp('#############################')
disp('L_beta = ')
disp(num2str(L_beta_exp))
disp('#############################')
disp('Lhat_beta = ')
disp(num2str(Lhat_beta_exp))
%%%%%%%%%%%%%%%%%%% RATIOS
disp('############# RATIOS ###############################')
disp('############# Galperin = 2 ############')
disp('Lhat_beta/LM = ')
disp(num2str(Lhat_beta_exp./LM_exp))
disp('############# Galperin = 4.6 ############')
disp('LM/L_beta = ')
disp(num2str(LM_exp./L_beta_exp))
%%%%%%%%%%%%%%%%%%% n intersection
disp('############# nth zeros ############')
disp('nth = ')
disp(num2str(n_intersec))
