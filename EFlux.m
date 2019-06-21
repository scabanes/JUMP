%%
%%% Script per il calcolo del flusso di energia
%%-------------------------------------------------------------------------
%__________________________________________________________________________
%% INPUTS:
% - U, V and Vort: are the velocity dans vorticity field.
%   [Fields] = (nodi*nodi, time)
% - nodi: is the length of the squared cartesian grid.
% - Vl: is a vector with the scale in grid point dimension you want to
%   filter on. It is of size N_l. The associated physical scale are
%   converted in cm using gridtocm.
% - passo_fr: we can choose the time step.
%__________________________________________________________________________
%% OUTPUTS:
% - FlussoEnergia_time_l & FlussoEnstrofia_time_l
%   This are the fluxes on the whole cartesian grid, for all times called 
%   and for the request scales in Vl. Scales are converted in cntimeters
%   using gridtocm and in wavenumbers using k = 2pi/(l*gridtocm).
%__________________________________________________________________________
clear s H Hc% give in Create_Grid_pol_Vasca and can be used in the following.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                   LOAD:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[U,b]=loadmtx([roots,Name,NameU]);
[V,b]=loadmtx([roots,Name,NameV]);
[Vort,b]=loadmtx([roots,Name,NameVort]);
% We work on a squared cartesian grid. nodi is the length of it. 
nodi = b(1)^0.5;
passo_fr = 5;
gridtocm = (data7*cmpx)/nodi;
Vl = [0.4:0.4:4,5:nodi];
N_l = length(Vl);

%  ###########################################################################################################################
%  #################################################### Filtraggio ###########################################################
%  ###########################################################################################################################
% [Field_fil_time_l,t_max] = GaussianFilter(Field,Vl,nodi,Tmax,Itime,passo_fr)
[U_fil_time_l,t_max] = GaussianFilter(U,Vl,N_l,nodi,Tmax,Itime,passo_fr);
[V_fil_time_l,t_max] = GaussianFilter(V,Vl,N_l,nodi,Tmax,Itime,passo_fr);
[UU_fil_time_l,t_max] = GaussianFilter(U.*U,Vl,N_l,nodi,Tmax,Itime,passo_fr);
[VV_fil_time_l,t_max] = GaussianFilter(V.*V,Vl,N_l,nodi,Tmax,Itime,passo_fr);
[UV_fil_time_l,t_max] = GaussianFilter(U.*V,Vl,N_l,nodi,Tmax,Itime,passo_fr);
[Vort_fil_time_l,t_max] = GaussianFilter(Vort,Vl,N_l,nodi,Tmax,Itime,passo_fr);
[VortU_fil_time_l,t_max] = GaussianFilter(Vort.*U,Vl,N_l,nodi,Tmax,Itime,passo_fr);
[VortV_fil_time_l,t_max] = GaussianFilter(Vort.*V,Vl,N_l,nodi,Tmax,Itime,passo_fr);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                      INITIATE MATRICES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Du_Dx  = zeros (nodi,nodi);
Dv_Dx  = zeros (nodi,nodi);
Du_Dy = zeros (nodi,nodi);
Dv_Dy = zeros (nodi,nodi);
%
DVort_Dx  = zeros (nodi,nodi);
DVort_Dy = zeros (nodi,nodi);
%
Flusso_Energia = zeros (nodi,nodi);
FlussoEnergia_time_l = zeros (nodi*nodi,t_max,N_l);
Flusso_Enstrofia = zeros (nodi,nodi);
FlussoEnstrofia_time_l = zeros (nodi*nodi,t_max,N_l);
% dx is dimensional, here in cm dans corresponds to the length of a grid
% point.
dx = gridtocm;
for t = 1:t_max
    
    disp(t)
    
    for  il = 1 : N_l
    
        % We reshape for the energy flux      
        U_fil_i_l = reshape (U_fil_time_l(:,t,il),nodi,nodi);
        V_fil_i_l = reshape (V_fil_time_l(:,t,il),nodi,nodi);
        UU_fil_i_l = reshape (UU_fil_time_l(:,t,il),nodi,nodi);
        VV_fil_i_l = reshape (VV_fil_time_l(:,t,il),nodi,nodi);
        UV_fil_i_l = reshape (UV_fil_time_l(:,t,il),nodi,nodi);
        % We reshape for the enstrophy flux
        Vort_fil_i_l = reshape (Vort_fil_time_l(:,t,il),nodi,nodi);
        VortV_fil_i_l = reshape (VortV_fil_time_l(:,t,il),nodi,nodi);
        VortU_fil_i_l = reshape (VortU_fil_time_l(:,t,il),nodi,nodi);
        
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                          DERIVATE D/Dx:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% ciclo per la derivata rispetto ad x
        for q = 1:nodi-1
            
            if q == 1
                
                Du_Dx (:,q) = (U_fil_i_l(:,q+1)-U_fil_i_l(:,q))/dx;
                Dv_Dx (:,q) = (V_fil_i_l(:,q+1)-V_fil_i_l(:,q))/dx;
                %
                DVort_Dx (:,q) = (Vort_fil_i_l(:,q+1)-Vort_fil_i_l(:,q))/dx;

               
            else
                
            Du_Dx (:,q) = (U_fil_i_l(:,q+1)-U_fil_i_l(:,q-1))/(2*dx);
            Dv_Dx (:,q) = (V_fil_i_l(:,q+1)-V_fil_i_l(:,q-1))/(2*dx);
            %
            DVort_Dx (:,q) = (Vort_fil_i_l(:,q+1)-Vort_fil_i_l(:,q-1))/(2*dx);

            
            end
            
        end
        
        for q = nodi
            
              Du_Dx (:,q) = (U_fil_i_l(:,q)-U_fil_i_l(:,q-1))/dx;
              Dv_Dx (:,q) = (V_fil_i_l(:,q)-V_fil_i_l(:,q-1))/dx;
              %
              DVort_Dx (:,q) = (Vort_fil_i_l(:,q)-Vort_fil_i_l(:,q-1))/dx;
            
        end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                          DERIVATE D/Dy:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% ciclo per la derivata rispetto a y..dx=dy
        for p = 1: nodi-1
            
            if p == 1
                
                Du_Dy (p,:) = (U_fil_i_l (p+1,:)-U_fil_i_l (p,:))/dx;
                Dv_Dy (p,:) = (V_fil_i_l (p+1,:)-V_fil_i_l (p,:))/dx;
                %
                DVort_Dy (p,:) = (Vort_fil_i_l (p+1,:)-Vort_fil_i_l (p,:))/dx;
               
            else
                
            Du_Dy (p,:) = (U_fil_i_l (p+1,:)-U_fil_i_l (p-1,:))/(2*dx);
            Dv_Dy (p,:) = (V_fil_i_l (p+1,:)-V_fil_i_l (p-1,:))/(2*dx);
            %
            DVort_Dy (p,:) = (Vort_fil_i_l (p+1,:)-Vort_fil_i_l (p-1,:))/(2*dx);
            
            end
            
        end
        
        for p = nodi
            
              Du_Dy (p,:) = (U_fil_i_l (p,:)-U_fil_i_l (p-1,:))/dx;
              Dv_Dy (p,:) = (V_fil_i_l (p,:)-V_fil_i_l (p-1,:))/dx;
              %
              DVort_Dy (p,:) = (Vort_fil_i_l (p,:)-Vort_fil_i_l (p-1,:))/dx;
            
        end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                 FLUSSI:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%:.............................................................
        %%% Inicio ciclo per calcolare il flusso
        
        
        for i = 1: nodi
            for j = 1:nodi
                
                Flusso_Energia (i,j) = -(Du_Dx(i,j)*(UU_fil_i_l(i,j)-(U_fil_i_l(i,j)*U_fil_i_l(i,j)))+...
                                         Dv_Dy(i,j)*(VV_fil_i_l(i,j)-(V_fil_i_l(i,j)*V_fil_i_l(i,j)))+...
                                        (Du_Dy(i,j)+Dv_Dx(i,j))*(UV_fil_i_l(i,j)-(U_fil_i_l(i,j)*V_fil_i_l(i,j))));
                
                Flusso_Enstrofia (i,j) = -(DVort_Dx(i,j)*(VortU_fil_i_l(i,j)-(U_fil_i_l(i,j)*Vort_fil_i_l(i,j)))+... 
                                           DVort_Dy(i,j)*(VortV_fil_i_l(i,j)-(V_fil_i_l(i,j)*Vort_fil_i_l(i,j)))); 
                             
            end
        end
        
        FlussoEnergia_time_l (:,t,il) = Flusso_Energia(:);
        FlussoEnstrofia_time_l (:,t,il) = Flusso_Enstrofia(:);
        
    end
end
FF = FlussoEnergia_time_l(:);
HH = FlussoEnstrofia_time_l(:);
% figure; hold on
% plot(Vl.*gridtocm, squeeze(mean(FlussoEnergia_time_l)))
% plot(Vl.*gridtocm, squeeze(mean(FlussoEnstrofia_time_l)))

%% ###########################################################################################################################
%  ###########################################################################################################################
%                                                   SAVE
%  ###########################################################################################################################
%  ###########################################################################################################################
fileout = [roots,Name,'/Energy_flux'] ;
filename = sprintf('%s.mtx',fileout);
 
fid = fopen(filename,'wb');
fwrite(fid,size(FF,1),'ulong');
fwrite(fid,size(FF,2),'ulong');
fwrite(fid,FF(:),'float');
fclose(fid);
%
fileout = [roots,Name,'/Enstrophy_flux'] ;
filename = sprintf('%s.mtx',fileout);
 
fid = fopen(filename,'wb');
fwrite(fid,size(HH,1),'ulong');
fwrite(fid,size(HH,2),'ulong');
fwrite(fid,HH(:),'float');
fclose(fid);
 