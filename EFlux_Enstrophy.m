%%% Script per il flusso di energia
clear s H Hc% give in Create_Grid_pol_Vasca and can be used in the following.

EFlux_GaussFilter_U_V

EFlux_GaussFilter_Vort_UVort_VVort

%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                  Derivate & Flussi
%  ###########################################################################################################################
%  ###########################################################################################################################
% 
% nodi = 128;
% 
% r_max = 18;
% 
% L_max = 2*r_max;
% N_l = L_max/dl;
% 
% totaltime = b(1)/(nodi*nodi)/N_l;
% 
% 
% 
% %%reashapre tutti i campi
% 
% Vort_fil = reshape (Vort_fil,nodi*nodi,totaltime,N_l);
% VortU_fil = reshape (VortU_fil,nodi*nodi,totaltime,N_l);
% VortV_fil = reshape (VortV_fil,nodi*nodi,totaltime,N_l);
% U_fil = reshape (U_fil,nodi*nodi,totaltime,N_l);
% V_fil = reshape (V_fil,nodi*nodi,totaltime,N_l);
% 
dx=(2*r_max_cm)/nodi;
Tmax =t_max;

% %%inizializzo le derivate e il flusso
 DVort_Dx  = zeros (nodi,nodi);
 DVort_Dy = zeros (nodi,nodi);
 
 
 Flusso = zeros (nodi,nodi);
 Flusso_time_l = zeros (nodi*nodi,Tmax,N_l);
 
for t = 1:Tmax
    
    disp(t)
    
     for  s = 1 : N_l
    
    
        Vort_fil_i_l = reshape (Vort_fil_time_l(:,t,s),nodi,nodi);
        VortV_fil_i_l = reshape (VortV_fil_time_l(:,t,s),nodi,nodi);
        VortU_fil_i_l = reshape (VortU_fil_time_l(:,t,s),nodi,nodi);
        U_fil_i_l = reshape (U_fil_time_l(:,t,s),nodi,nodi);
        V_fil_i_l = reshape (V_fil_time_l(:,t,s),nodi,nodi);       
        
        %%%%% ciclo per la derivata rispetto ad x
        for q = 1:nodi-1
            
            if q == 1
                
                DVort_Dx (:,q) = (Vort_fil_i_l(:,q+1)-Vort_fil_i_l(:,q))/dx;
                               
            else
                
            DVort_Dx (:,q) = (Vort_fil_i_l(:,q+1)-Vort_fil_i_l(:,q-1))/(2*dx);
           
            
            end
            
        end
        
        for q = nodi
            
              DVort_Dx (:,q) = (Vort_fil_i_l(:,q)-Vort_fil_i_l(:,q-1))/dx;
             
            
        end
        
        %%%%%%% ciclo per la derivata rispetto a y..dx=dy
        for p = 1: nodi-1
            
            if p == 1
                
                DVort_Dy (p,:) = (Vort_fil_i_l (p+1,:)-Vort_fil_i_l (p,:))/dx;
                              
            else
                
            DVort_Dy (p,:) = (Vort_fil_i_l (p+1,:)-Vort_fil_i_l (p-1,:))/(2*dx);
           
            end
            
        end
        
        for p = nodi
            
              DVort_Dy (p,:) = (Vort_fil_i_l (p,:)-Vort_fil_i_l (p-1,:))/dx;
             
        end
        
        %%%%:.............................................................
        %%% Inicio ciclo per calcolare il flusso
        
        
        for i = 1: nodi
            for j = 1:nodi
                
%                 Flusso (i,j) = -(DVort_Dx(i,j)*(VortU_fil_i_l(i,j)+VortV_fil_i_l(i,j)-(2*U_fil_i_l(i,j)*Vort_fil_i_l(i,j)))+... 
%                                 DVort_Dy(i,j)*(VortU_fil_i_l(i,j)+VortV_fil_i_l(i,j)-(2*V_fil_i_l(i,j)*Vort_fil_i_l(i,j))));
                Flusso (i,j) = -(DVort_Dx(i,j)*(VortU_fil_i_l(i,j)-(U_fil_i_l(i,j)*Vort_fil_i_l(i,j)))+... 
                                DVort_Dy(i,j)*(VortV_fil_i_l(i,j)-(V_fil_i_l(i,j)*Vort_fil_i_l(i,j)))); 
                            
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
fileout = [roots,Name,'/Enstrophy_flux_dl_' num2str(dl)] ;
filename = sprintf('%s.mtx',fileout);
 
fid = fopen(filename,'wb');
fwrite(fid,size(FF,1),'ulong');
fwrite(fid,size(FF,2),'ulong');
fwrite(fid,FF(:),'float');
fclose(fid);
 

