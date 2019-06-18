%%% Costruzione del filtro di Gaussiano e filtraggio dei campi di velocit� 
%%% in base alla scala di cut-off.
[U,b]=loadmtx([roots,Name,NameU]);
[V,b]=loadmtx([roots,Name,NameV]);
% We are on a squared cartesian grid with a size of:
nodi = b(1)^0.5;
%-->> Maps
figure; hold on
title('initial field')
contourf(reshape(U(:,1),nodi,nodi),40,'LineStyle','none')
colorbar
%%% considero la diagonale del quadrato
L_max = 2*r_max_cm;
%N_l = L_max/dl;
N_l = round(4/dl + (L_max-4));

%%passo temporale per l'analisi
passo_fr = 5;
t_max= floor(Tmax/passo_fr);

s=0;

%%% U matrice 3D dei campi di u filtrati le colonne sono i tempi, le righe
%%% lo spazio 128*128 srotolato, la terza dimensione sono le scale di cut-off

U_fil_time_l = zeros (nodi*nodi,t_max,N_l);
V_fil_time_l = zeros (nodi*nodi,t_max,N_l);


for l = dl : dl:4

 
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

   % disp(t)
  %vorticit� normalizzata per f0
    U_i =reshape(U(:,t),nodi,nodi);
    V_i =reshape(V(:,t),nodi,nodi);
     
   % disp(t)
  
   U_fil = filter2(h,U_i);
   V_fil = filter2(h,V_i);

   U_fil_time_l (:,appo,s) = U_fil (:);
   V_fil_time_l (:,appo,s) = V_fil (:);
   
   

end

end

for l = 5: L_max

 
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

   % disp(t)
  
    U_i =reshape(U(:,t),nodi,nodi);
    V_i =reshape(V(:,t),nodi,nodi);
     
   % disp(t)
  
   U_fil = filter2(h,U_i);
   V_fil = filter2(h,V_i);

   U_fil_time_l (:,appo,s) = U_fil (:);
   V_fil_time_l (:,appo,s) = V_fil (:);
   
   

end

end

figure; hold on
title('filtred field')
contourf(U_fil,40,'LineStyle','none')
colorbar

%%% decommentare se si fa girare da solo o se si vogliono salvare i campi
%%% filtrati di U e V

%  U_FIL = U_fil_time_l(:);
%  
%  V_FIL = V_fil_time_l(:);
%  
%  fileout = ['U' num2str(n) '_fil_dl_1'] ;
%  
%  filename = sprintf('%s.mtx',fileout);
%  
%  fid = fopen(filename,'wb');
%  fwrite(fid,size(U_FIL,1),'ulong');
%  fwrite(fid,size(U_FIL,2),'ulong');
%  fwrite(fid,U_FIL(:),'float');
%  fclose(fid);
% 
%  
% fileout = ['V' num2str(n) '_fil_dl_1'] ;
%  
%  filename = sprintf('%s.mtx',fileout);
%  
%  fid = fopen(filename,'wb');
%  fwrite(fid,size(V_FIL,1),'ulong');
%  fwrite(fid,size(V_FIL,2),'ulong');
%  fwrite(fid,V_FIL(:),'float');
%  fclose(fid);
%  