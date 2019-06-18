%% Script che lavora a seguito dello script Filtraggio_Gauss_U_V_OK.m per
%% il calcolo dei campi filtrati di U*U, V*V, U*V
%% ---------------------------------------------------------------------

UV = U.*V;

s=0;

UU_fil_time_l =  zeros (nodi*nodi,t_max,N_l);
VV_fil_time_l =  zeros (nodi*nodi,t_max,N_l);
UV_fil_time_l =  zeros (nodi*nodi,t_max,N_l);


for l = dl:dl:4

    disp (l)
    s = s+1;
    
    % calcolo la deviazione standard
    sig = l/(12)^0.5;
   
      
  if sig <= 1
    
    h = fspecial ('gaussian',5,sig);
   
  else
    
    dim_max = round(2*(3*sig));
    
    if mod(dim_max,2)  ==  0
  
       %disp('dimensione pari');
   
         W = dim_max + 1;
        
    h = fspecial ('gaussian',W,sig);
    
    else
        
     W = dim_max  ;
     h = fspecial ('gaussian',W,sig);
   
    end
    
     
  end

%end

appo=0;

for t = Itime : passo_fr : Itime+Tmax-1
    
    appo=appo+1;

   % disp(t)
  
  U_i =reshape(U(:,t),nodi,nodi);
  U_i = U_i.^2;
  
  V_i =reshape(V(:,t),nodi,nodi);
  V_i = V_i.^2;
  
  UV_i = reshape(UV(:,t),nodi,nodi);
    
  UU_fil = filter2(h,U_i);
  VV_fil = filter2(h,V_i);
  UV_fil = filter2(h,UV_i);


   UU_fil_time_l (:,appo,s) = UU_fil (:);
   VV_fil_time_l (:,appo,s) = VV_fil (:);
   UV_fil_time_l (:,appo,s) = UV_fil (:);  

end

end

for l = 5:L_max

    disp (l)
    s = s+1;
    
    % calcolo la deviazione standard
    sig = l/(12)^0.5;
   
      
  if sig <= 1
    
    h = fspecial ('gaussian',5,sig);
   
  else
    
    dim_max = round(2*(3*sig));
    
    if mod(dim_max,2)  ==  0
  
       %disp('dimensione pari');
   
         W = dim_max + 1;
        
    h = fspecial ('gaussian',W,sig);
    
    else
        
     W = dim_max  ;
     h = fspecial ('gaussian',W,sig);
   
    end
    
     
  end

%end

appo=0;

for t = Itime : passo_fr : Itime+Tmax-1
    
    appo=appo+1;

   % disp(t)
  %vorticit� normalizzata per f0
  U_i = reshape(U(:,t),nodi,nodi);
  U_i = U_i.^2;
  
  V_i = reshape(V(:,t),nodi,nodi);
  V_i = V_i.^2;
  
  UV_i = reshape(UV(:,t),nodi,nodi);
  
  
  UU_fil = filter2(h,U_i);
  VV_fil = filter2(h,V_i);
  UV_fil = filter2(h,UV_i);


   UU_fil_time_l (:,appo,s) = UU_fil (:);
   VV_fil_time_l (:,appo,s) = VV_fil (:);
   UV_fil_time_l (:,appo,s) = UV_fil (:);  

end

end

%  U_FIL = U_fil_time_l(:);
%  
%  V_FIL = V_fil_time_l(:);
%  
%   UV_FIL = UV_fil_time_l(:);
%   
%  
%  fileout = ['UU' num2str(n) '_fi_dl_1'] ;
%  
%  filename = sprintf('%s.mtx',fileout);
%  fid = fopen(filename,'wb');
%  fwrite(fid,size(U_FIL,1),'ulong');
%  fwrite(fid,size(U_FIL,2),'ulong');
%  fwrite(fid,U_FIL(:),'float');
%  fclose(fid);
% 
%  
%  fileout = ['VV' num2str(n) '_fi_dl_1'] ;
%  
%  filename = sprintf('%s.mtx',fileout);
%  
%  fid = fopen(filename,'wb');
%  fwrite(fid,size(V_FIL,1),'ulong');
%  fwrite(fid,size(V_FIL,2),'ulong');
%  fwrite(fid,V_FIL(:),'float');
%  fclose(fid);
%  
%  fileout = ['UV' num2str(n) '_fi_dl_1'] ;
%  
%  filename = sprintf('%s.mtx',fileout);
%  
%  fid = fopen(filename,'wb');
%  fwrite(fid,size(UV_FIL,1),'ulong');
%  fwrite(fid,size(UV_FIL,2),'ulong');
%  fwrite(fid,UV_FIL(:),'float');
%  fclose(fid);
%  