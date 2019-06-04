%%%%% Programma per interpolare il dato sparso contenuto nel FES su griglia
%%%%% polare creata però in px con lo script Create_Grid_pol_OK.m
%%%%% PER ESPERIMENTI4 e %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Input:   - il fes da interpolare
%% Output:  - U_pol, V_pol, componenti della velocità (u e v cartesiane) sui punti di una griglia polare       

%% Nota: - le velocità sono in cm/s
%        - sono salvate in questo modo: ogni colonna è un istante termporale 
%%         da reshepare così U_t = reshape (U_pol(:,t),Nraggi,Ncerchi), 
%%         U_t sarà una matrice in cui ogni colonna rappresenta un cerchio a distanza r dal centro
%%         le righe rappresentano le distanze angolari
close all
clear all
%%-------------------------------------------------------------------------
%% creo la griglia polare in px definendo il numero di raggi e di cerchi.
%% La griglia andrà da 0 a Roiwidth (es. 900 px)
%**********************************************************
%           PARAMETRI DA SCEGLIERE
%**********************************************************
cexp = 2;
expt = {'5';'6';'7';'8';'9';'3'}
roots = '../Dati-Torino/';
output = 'outputs'
ifsave = 1; %if 1 files are saved
%**********************************************************
% ##################################################################################################################################################
%                                                                                                                                 DATI INFORMAZIONI:
% ##################################################################################################################################################
EXPT_infos
% ##################################################################################################################################################
%                                                                                                                                             GRID :
% ##################################################################################################################################################
Create_Grid_pol_Torino_m
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nFrames = ((endi{cexp,2}-starti{cexp,2})/step + 1);
% nFrames=2150;  %%%solo exp 8 che ha un frame mancante!!!!!
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% inizializzazione
Vz_t = zeros(Nraggi*Ncerchi,nFrames);
Vr_t = zeros(Nraggi*Ncerchi,nFrames);
U_int_ii_t = zeros(Nraggi*Ncerchi,nFrames);
V_int_ii_t = zeros(Nraggi*Ncerchi,nFrames);
%-------------------------------------------------------------------------
% ##################################################################################################################################################
%                                                                                                                                   INTERPOLAZIONE :
% ##################################################################################################################################################
for t=1:nFrames
  %  lenght=fread(fid2,1,'int32');
  t
  
imagei(t) = starti{cexp,2} + step * (t-1);
imageo(t) = starto{cexp,2} + step * (t-1);

% file_inner = "%s/%s/inner_%s.sback.civ.calib/fig_%d-%d.nc" % (prefix,names[r],ftype[r],imagei[t],imagei[t]+gapi[r])
file_inner = [roots,names{cexp,2},'/inner_',ftype{cexp,2},'.sback.civ.calib/fig_',num2str(imagei(t)),'-',num2str(imagei(t)+gapi{cexp,2}),'.nc'];
file_outer = [roots,names{cexp,2},'/outer_',ftype{cexp,2},'.sback.civ.calib/fig_',num2str(imageo(t)),'-',num2str(imageo(t)+gapo{cexp,2}),'.nc'];

% Load all files
U_out = double(ncread(file_outer,'Civ2_U'));
V_out = double(ncread(file_outer,'Civ2_V'));
X_out = double(ncread(file_outer,'Civ2_X'));
Y_out = double(ncread(file_outer,'Civ2_Y'));
U_in = double(ncread(file_inner,'Civ2_U'));
V_in = double(ncread(file_inner,'Civ2_V'));
X_in = double(ncread(file_inner,'Civ2_X'));
Y_in = double(ncread(file_inner,'Civ2_Y'));
% concatenate in one vector inner and outer camera
U = [U_out;U_in];
V = [V_out;V_in];
X = [X_out;X_in];
Y = [Y_out;Y_in];
L = length(X);
%quiver (X,Y,U, V);
%----------------------------------------------
% ==> Toward interpolation:
% reshape the query polar grid in a vector
SS=size(Grid_Xp_m);   %%verifica questo!
VGrid_Xp_m = reshape(Grid_Xp_m,[SS(1)*SS(2),1]);
VGrid_Yp_m = reshape(Grid_Yp_m,[SS(1)*SS(2),1]);
% Find the polygonal boundary of the data defined
% at coordinates X and Y. Note indice 1 is the 
% qccurqncy of the boundary
k = boundary(X,Y,1);
% Find in vectors of query polar grid the coordinates
% that are inside the bordered polygon -> logical vector in
in = inpolygon(VGrid_Xp_m,VGrid_Yp_m,X(k),Y(k));
% Interpolation on the query polqr grid in the polygonal
% data boundary.

U_int_i = griddata (X,Y,U,VGrid_Xp_m(in),VGrid_Yp_m(in),'cubic');
V_int_i= griddata (X,Y,V,VGrid_Xp_m(in),VGrid_Yp_m(in),'cubic');

% From vector to a Nraggi*Nerchi array
MU_z = zeros(size(Grid_Xp_m));
MU_r = zeros(size(Grid_Xp_m));

ii=1;
jj=1;
for i = 1:SS(2)
    for j = 1:SS(1)
        if(in(jj) == 1)
            MU_z(j,i) = ((-U_int_i(ii))*sin(Tpoints(j))+(V_int_i(ii))*cos(Tpoints(j)));
            MU_r(j,i) =  U_int_i(ii)*cos(Tpoints(j))+(V_int_i(ii))*sin(Tpoints(j));
            ii = ii+1;
        else
            MU_z(j,i) = nan;
            MU_r(j,i) = nan;
        end
        jj = jj+1;
    end
end
Vz_t (:,t) = MU_z(:);
Vr_t (:,t) = MU_r(:);

%per scrivere le componenti cartesiane su griglia polare%%%%%
U_int_ii=zeros(size(Grid_Xp_m));
V_int_ii=zeros(size(Grid_Xp_m));  

ii=1;
jj=1;
for i = 1:SS(2)
    for j = 1:SS(1)
        if(in(jj) == 1)
            U_int_ii(j,i) = U_int_i(ii);
            V_int_ii(j,i) = V_int_i(ii);
            ii = ii+1;
        else
            U_int_ii(j,i) = nan;
            V_int_ii(j,i) = nan;
        end
        jj = jj+1;
    end
end
U_int_ii_t(:,t) = U_int_ii(:);
V_int_ii_t(:,t) = V_int_ii(:);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% size(X);
% SS=size(Grid_Xp_m);
% VGrid_Xp_m = reshape(Grid_Xp_m,[SS(1)*SS(2),1]);
% VGrid_Yp_m = reshape(Grid_Yp_m,[SS(1)*SS(2),1]);
% k = boundary(X,Y,1);
% in = inpolygon(VGrid_Xp_m,VGrid_Yp_m,X(k),Y(k));
% 
% Grid_Xp_m(in)
% % %     figure
% % %     quiver(X,Y,U,V,'k')
% % %     hold on
% % %  %plot(X(k),Y(k),'k')
% % %     pcolor(Grid_Xp_m,Grid_Yp_m,MU_z); shading interp
% % %  %quiver(VGrid_Xp_m(in),VGrid_Yp_m(in),ones(size(VGrid_Yp_m(in))),ones(size(VGrid_Yp_m(in))),'r')
% 
%% dimensionalizzazione delle velocità in cm/s
%U_int = U_int *F *cmpx;
%V_int = V_int *F *cmpx;

%% salvataggio delle componenti di U e V su griglia polare in .mtx
% ##################################################################################################################################################
%                                                                                                                                            PLOTS :
% ##################################################################################################################################################
%quiver(Grid_Xp_m,Grid_Yp_m,U_int_i, V_int_i)
 figure
 contourf(Grid_Xp_m,Grid_Yp_m,U_int_ii, 100, ':');
 hold on 
% quiver(Grid_Xp_m(idx),Grid_Yp_m(idx),U_int_i,V_int_i,'r')
% hold on
quiver(X,Y,U,V,'r')
grid on
% ##################################################################################################################################################
%                                                                                                                                    SALVIAMO DATI :
% ##################################################################################################################################################
% save([roots,names{cexp,2},'/',output,'/infos.mat'],'Nraggi','Ncerchi','dphi','dr','r','MaxRadius','Rpoints','Tpoints')
if (ifsave==1)
 if ~exist([roots,names{cexp,2},'/',output], 'dir')
       mkdir([roots,names{cexp,2},'/',output])
 end
%save polar component in .mtx files
%zonal
%fileout = [roots,'Vz09_540_180_2160'];
fileout = [roots,names{cexp,2},'/',output,'/Vz_',num2str(Nraggi),'_',num2str(Ncerchi),'_',num2str(nFrames)];
filename = sprintf('%s.mtx',fileout);
fid = fopen(filename,'wb');
fwrite(fid,size(Vz_t,1),'ulong');
fwrite(fid,size(Vz_t,2),'ulong');
fwrite(fid,Vz_t(:),'float');
fclose(fid);

% % % %radial
% % % 
fileout2 = [roots,names{cexp,2},'/',output,'/Vr_',num2str(Nraggi),'_',num2str(Ncerchi),'_',num2str(nFrames)];
filename2 = sprintf('%s.mtx',fileout2);
fid = fopen(filename2,'wb');
fwrite(fid,size(Vr_t,1),'ulong');
fwrite(fid,size(Vr_t,2),'ulong');
fwrite(fid,Vr_t(:),'float');
fclose(fid);
% % % 
% % % U cartesiana su griglia polare
% fileout3 = ['U9_pol_540_180_2160' ];
fileout3 = [roots,names{cexp,2},'/',output,'/U_pol_',num2str(Nraggi),'_',num2str(Ncerchi),'_',num2str(nFrames)];
filename3 = sprintf('%s.mtx',fileout3);
fid = fopen(filename3,'wb');
fwrite(fid,size(U_int_ii_t,1),'ulong');
fwrite(fid,size(U_int_ii_t,2),'ulong');
fwrite(fid,U_int_ii_t(:),'float');
fclose(fid);
% % % 
% % %  V cartesiana su griglia polare
% fileout4 = ['V9_pol_540_180_2160' ];
fileout4 = [roots,names{cexp,2},'/',output,'/V_pol_',num2str(Nraggi),'_',num2str(Ncerchi),'_',num2str(nFrames)];
filename4 = sprintf('%s.mtx',fileout4);
fid = fopen(filename4,'wb');
fwrite(fid,size(V_int_ii_t,1),'ulong');
fwrite(fid,size(V_int_ii_t,2),'ulong');
fwrite(fid,V_int_ii_t,'float');
fclose(fid);
% % % 
end
save([roots,names{cexp,2},'/',output,'/infos.mat'],'Nraggi','Ncerchi','dphi','dr','r','MaxRadius','Rpoints','Tpoints')