%%
%%% Script per il calcolo del flusso di energia
%%-------------------------------------------------------------------------
%% OUTPUT: un vettore srotolato con i valori del flusso di energia in ogni punto dello spazio, 
%          per ogni istante temporale (analisi di 1 fr ogni 5), per ogni
%          scala di cut-off.
%          da reshapare cos�
%          F_Energy = reshape(F_Energy,nodi*nodi,round(totaltime),N_l);
%          dove nodi = 128
%          totaltime lo ricavo nelloscript dei plot
%          N_l � il numero massimo di scale considerate per il calcolo del
%          flusso
%%% servono nella stessa cartella i primi due script mensionati. 
clear s H Hc% give in Create_Grid_pol_Vasca and can be used in the following.

EFlux_GaussFilter_U_V

EFlux_GaussFilter_UU_VV_UV

%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                  Derivate & Flussi
%  ###########################################################################################################################
%  ###########################################################################################################################
% 
% %%reashapre tutti i campi
% 
% U_fil = reshape (U_fil,nodi*nodi,totaltime,N_l);
% UU_fil = reshape (UU_fil,nodi*nodi,totaltime,N_l);
% V_fil = reshape (V_fil,nodi*nodi,totaltime,N_l);
% VV_fil = reshape (VV_fil,nodi*nodi,totaltime,N_l);
% UV_fil = reshape (UV_fil,nodi*nodi,totaltime,N_l);

%%%------------------------------------------------------------------------
Tmax = t_max;
dx=(2*r_max_cm)/nodi;

%%inizializzo le derivate e il flusso
Du_Dx  = zeros (nodi,nodi);
Dv_Dx  = zeros (nodi,nodi);
Du_Dy = zeros (nodi,nodi);
Dv_Dy = zeros (nodi,nodi);

Flusso = zeros (nodi,nodi);
Flusso_time_l = zeros (nodi*nodi,Tmax,N_l);

s=0;

for t = 1:Tmax
    
    disp(t)
    
    for  s = 1 : N_l
    
               
        U_fil_i_l = reshape (U_fil_time_l(:,t,s),nodi,nodi);
        V_fil_i_l = reshape (V_fil_time_l(:,t,s),nodi,nodi);
        UU_fil_i_l = reshape (UU_fil_time_l(:,t,s),nodi,nodi);
        VV_fil_i_l = reshape (VV_fil_time_l(:,t,s),nodi,nodi);
        UV_fil_i_l = reshape (UV_fil_time_l(:,t,s),nodi,nodi);
        
        
        %%%%% ciclo per la derivata rispetto ad x
        for q = 1:nodi-1
            
            if q == 1
                
                Du_Dx (:,q) = (U_fil_i_l(:,q+1)-U_fil_i_l(:,q))/dx;
                Dv_Dx (:,q) = (V_fil_i_l(:,q+1)-V_fil_i_l(:,q))/dx;
               
            else
                
            Du_Dx (:,q) = (U_fil_i_l(:,q+1)-U_fil_i_l(:,q-1))/(2*dx);
            Dv_Dx (:,q) = (V_fil_i_l(:,q+1)-V_fil_i_l(:,q-1))/(2*dx);
            
            end
            
        end
        
        for q = nodi
            
              Du_Dx (:,q) = (U_fil_i_l(:,q)-U_fil_i_l(:,q-1))/dx;
              Dv_Dx (:,q) = (V_fil_i_l(:,q)-V_fil_i_l(:,q-1))/dx;
            
        end
        
        %%%%%%% ciclo per la derivata rispetto a y..dx=dy
        for p = 1: nodi-1
            
            if p == 1
                
                Du_Dy (p,:) = (U_fil_i_l (p+1,:)-U_fil_i_l (p,:))/dx;
                Dv_Dy (p,:) = (V_fil_i_l (p+1,:)-V_fil_i_l (p,:))/dx;
               
            else
                
            Du_Dy (p,:) = (U_fil_i_l (p+1,:)-U_fil_i_l (p-1,:))/(2*dx);
            Dv_Dy (p,:) = (V_fil_i_l (p+1,:)-V_fil_i_l (p-1,:))/(2*dx);
            
            end
            
        end
        
        for p = nodi
            
              Du_Dy (p,:) = (U_fil_i_l (p,:)-U_fil_i_l (p-1,:))/dx;
              Dv_Dy (p,:) = (V_fil_i_l (p,:)-V_fil_i_l (p-1,:))/dx;
            
        end
        
        %%%%:.............................................................
        %%% Inicio ciclo per calcolare il flusso
        
        
        for i = 1: nodi
            for j = 1:nodi
                
                Flusso (i,j) = -(Du_Dx(i,j)*(UU_fil_i_l(i,j)-(U_fil_i_l(i,j)*U_fil_i_l(i,j)))+...
                                 Dv_Dy(i,j)*(VV_fil_i_l(i,j)-(V_fil_i_l(i,j)*V_fil_i_l(i,j)))+...
                                 (Du_Dy(i,j)+Dv_Dx(i,j))*(UV_fil_i_l(i,j)-(U_fil_i_l(i,j)*V_fil_i_l(i,j))));
                             
            end
        end
        
        Flusso_time_l (:,t,s) = Flusso(:);
        
    end
end

FF = Flusso_time_l(:);

%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   SAVE
%  ###########################################################################################################################
%  ###########################################################################################################################
fileout = [roots,Name,'/Energy_flux_dl_' num2str(dl)] ;
filename = sprintf('%s.mtx',fileout);
 
fid = fopen(filename,'wb');
fwrite(fid,size(FF,1),'ulong');
fwrite(fid,size(FF,2),'ulong');
fwrite(fid,FF(:),'float');
fclose(fid);
 

