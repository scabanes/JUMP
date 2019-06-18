%%% Costruzione del filtro di Gaussiano e filtraggio dei campi di velocit�
%%% in base alla scala di cut-off..
%%% a partire da file mtx di velocit� dopo aver tagliato frame nel mtr

% clear all
[Vort,b]=loadmtx([roots,Name,NameVort]);

%%% U matrice 3D dei campi di u filtrati le colonne sono itempi, le righe
%%% lo spazio 128*128 srotolato, la terza dimensione sono le scale di cut-off

%%%costruzione del filtro
%%----------------------------------------------------------------------

s=0;


VortU = Vort.*U;
VortV = Vort.*V;

Vort_fil_time_l = zeros (nodi*nodi,t_max,N_l);
VortU_fil_time_l = zeros (nodi*nodi,t_max,N_l);
VortV_fil_time_l = zeros (nodi*nodi,t_max,N_l);


for l = dl : dl : 4 %L_max
 
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

for t = Itime : passo_fr : Itime+Tmax-1  % era commentato
    
    appo=appo+1;

    Vort_i =reshape(Vort(:,t),nodi,nodi);
    VortU_i =reshape(VortU(:,t),nodi,nodi);
    VortV_i =reshape(VortV(:,t),nodi,nodi);
     
   % disp(t)
  
   Vort_fil = filter2(h,Vort_i);
   VortU_fil = filter2(h,VortU_i);
   VortV_fil = filter2(h,VortV_i);

   Vort_fil_time_l(:,appo,s) = Vort_fil (:);
   VortU_fil_time_l (:,appo,s) = VortU_fil (:);
   VortV_fil_time_l (:,appo,s) = VortV_fil (:);
   
   

end

end

% seconda parte del ciclo
for l = 5 : L_max
 
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

for t = Itime : passo_fr : Itime+Tmax-1 % era commentato
    
    appo=appo+1;

    Vort_i =reshape(Vort(:,t),nodi,nodi);
    VortU_i =reshape(VortU(:,t),nodi,nodi);
    VortV_i =reshape(VortV(:,t),nodi,nodi);
     
   % disp(t)
  
   Vort_fil = filter2(h,Vort_i);
   VortU_fil = filter2(h,VortU_i);
   VortV_fil = filter2(h,VortV_i);

   Vort_fil_time_l(:,appo,s) = Vort_fil (:);
   VortU_fil_time_l (:,appo,s) = VortU_fil (:);
   VortV_fil_time_l (:,appo,s) = VortV_fil (:);
   
   

end

end

%  Vort_FIL = Vort_fil_time_l(:);
%  VortU_FIL = VortU_fil_time_l(:);
%  VortV_FIL = VortV_fil_time_l(:);
 
 %%% per reshapare U(nodi*nodi,t_max,l_max);
 
%  fileout = ['Vort' num2str(n) '_fil_dl_' num2str(dl) '_' num2str(I)] ;
%  
%  filename = sprintf('%s.mtx',fileout);
%  
%  fid = fopen(filename,'wb');
%  fwrite(fid,size(Vort_FIL ,1),'ulong');
%  fwrite(fid,size(Vort_FIL ,2),'ulong');
%  fwrite(fid,Vort_FIL (:),'float');
%  fclose(fid);
%  
%  fileout = ['VortU' num2str(n) '_fil_dl_' num2str(dl) '_' num2str(I)] ;
%  
%  filename = sprintf('%s.mtx',fileout);
%  
%  fid = fopen(filename,'wb');
%  fwrite(fid,size(VortU_FIL ,1),'ulong');
%  fwrite(fid,size(VortU_FIL ,2),'ulong');
%  fwrite(fid,VortU_FIL (:),'float');
%  fclose(fid);
% 
%  
% fileout = ['VortV' num2str(n) '_fil_dl_' num2str(dl) '_' num2str(I)] ;
%  
%  filename = sprintf('%s.mtx',fileout);
%  
%  fid = fopen(filename,'wb');
%  fwrite(fid,size(VortV_FIL ,1),'ulong');
%  fwrite(fid,size(VortV_FIL ,2),'ulong');
%  fwrite(fid,VortV_FIL (:),'float');
%  fclose(fid);
%  