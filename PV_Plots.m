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
[PV_t,b]=loadmtx([roots,'/VortPot_time_',num2str(Nti),'_',num2str(Nri),'_',num2str(nFrames)]);
[Vz_t,b] = loadmtx([roots,NameVt]);


figure; hold on
for it=1:1000
PV=reshape(PV_t(:,it),Nti,Nri);
Vz=reshape(Vz_t(:,it),Nti,Nri);

plot(PV(300,:),r)
plot(Vz(300,:),r,'r')
ylim([12 28])
end